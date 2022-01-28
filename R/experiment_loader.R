###############################################################################
# instantiate a sampler object for a given experiment
# note: the fields in "tuning_results.txt" files are
#   sd_thresh min_mass   sd_loglik time_per_loglik varmat_mult total_sampling_time ess_mean ess_per_min
###############################################################################

#######################################
# utilities common to all experiments
#######################################

# dispatcher, keep updated with experiments
#' @export
get_sampler = function(exp_name,reg_ts,gtp_solver,tuned_pars=TRUE,
                       use_theta_true=TRUE,...){
  caller = switch(
    EXPR = exp_name,
    GHS2017_LV  = get_sampler_GHS2017_LV,
    GHS2017_SIR = get_sampler_GHS2017_SIR,
    SG2019_LV20 = get_sampler_SG2019_LV,
    SG2019_Sch  = get_sampler_SG2019_Sch,
    MMc         = get_sampler_MMc,
    stop("Unknown experiment."))
  res_list = do.call(caller,list(reg_ts=reg_ts))
  alist=c(list(exp_name=exp_name,gtp_solver=gtp_solver),res_list$pars,list(...))
  if(tuned_pars){
    model_string = paste0(
      sprintf("%s_%s%s",exp_name,ifelse(reg_ts,"reg_ts_",""), gtp_solver))
    res_path = file.path("extdata","MH_sampler_tuning",model_string)
    alist$varmat = load_varmat( # read variance matrix
      path_to_varmat = system.file(res_path, "varmat.txt", package = "ctmc3"),
      par_names = alist$par_names)
    alist$theta_0=scan(
      system.file(res_path, "theta_0.txt", package = "ctmc3"),quiet=TRUE)
    names(alist$theta_0)=alist$par_names
    min_mass=scan(
      system.file(res_path, "tuning_results.txt", package = "ctmc3"),
      quiet=TRUE)[2L]
    sampler = do.call("new",alist,envir = res_list$sampler_constructor)
    sampler$study_convergence()
    sampler$set_AST(min_mass=min_mass)
  }else{
    alist$theta_0=if(use_theta_true) alist$theta_true else get_rnd_theta(alist$theta_true)
    sampler = do.call("new",alist,envir = res_list$sampler_constructor)
  }
  return(sampler)
}

# random initial point with the same scale as true parameter
get_rnd_theta=function(theta_true){mean(theta_true)*rexp(length(theta_true))}

# util for loading data into our format
load_exp_data = function(path_to_data){
  dta_raw = read.table(path_to_data,sep = " ")
  n_spec = ncol(dta_raw)-1L; n_obs = nrow(dta_raw)-1L
  dta = list(t = as.matrix(dta_raw[,1L]),
             s = as.matrix(unname(dta_raw[,-1L]),ncol=n_spec))
  dta$dt = as.matrix(c(NA_real_,diff(dta$t)))
  dta$s_pre = rbind(rep(NA_integer_, n_spec), dta$s[-(n_obs+1L),,drop=FALSE])
  dta = lapply(dta, function(v) v[-1L,,drop=FALSE])
  return(dta)
}

# util for loading variance matrices
load_varmat = function(path_to_varmat,par_names,path_to_alt){
  if(!file.exists(path_to_varmat)){
    cat("Tuned covariance unavailable, trying the alternative.\n")
    path_to_varmat = path_to_alt
  }
  varmat = as.matrix(read.table(file = path_to_varmat,sep = " "))
  dimnames(varmat) = replicate(2L,par_names,simplify = FALSE)
  return(varmat)
}

###############################################################################
# GHS2017 experiments
###############################################################################

#######################################
# Lotka-Volterra
#######################################

# prior logdensity used by the authors
GHS2017_LV_ldprior = function(theta){
  sum(dgamma(x=theta,4,10000,log=TRUE))
}

# extend ReactionNetwork to add LV propensities
ReactionNetworkLV = R6::R6Class(
  classname = "ReactionNetworkLV",
  inherit = ReactionNetwork,
  public = list(
    get_propensity = function(i,s){
      self$react_rates[i] * switch (i,
        prod(s), # birth of predator
        s[1L],   # death of predator
        s[2L],   # birth of prey
        prod(s)  # death of prey
      )
    }
  )
)
# extend to add specific prior
MHSamplerITSAdaptSTGHS2017LV = R6::R6Class(
  classname = "MHSamplerITSAdaptSTGHS2017LV",
  inherit = MHSamplerReactNetITSAdaptST,
  public = list(ldprior = GHS2017_LV_ldprior))
# extend to add specific prior
MHSamplerRTSAdaptSTGHS2017LV = R6::R6Class(
  classname = "MHSamplerRTSAdaptSTGHS2017LV",
  inherit = MHSamplerReactNetRTSAdaptST,
  public = list(ldprior = GHS2017_LV_ldprior))

# initialize sampler
get_sampler_GHS2017_LV = function(reg_ts=FALSE,...){
  # define CTMC model: Lotka-Volterra [predator(col1)-prey(col2)]
  updates = matrix(as.integer(c( 1, 0,
                                -1, 0,
                                 0, 1,
                                 0,-1)),
                   nrow=4L,byrow = TRUE)
  theta_true=c(0.0001, 0.0005, 0.0005, 0.0001) # trueVals in https://github.com/ageorgou/roulette/blob/master/experiments/predPreyMHRT.m
  par_names = c("Birth of predator","Death of predator",
                "Birth of prey","Death of prey")
  CTMC=ReactionNetworkLV$new( # instantiate a CTMC object
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "obsPredPreyRT_true", package = "ctmc3"))
  # dta = load_exp_data(file.path(".","data","obsPredPreyRT_true"))
  sampler_constructor = if(reg_ts) MHSamplerRTSAdaptSTGHS2017LV else MHSamplerITSAdaptSTGHS2017LV
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,
                 theta_true=theta_true)
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}

#######################################
# SIR
#######################################

# extend ReactionNetwork to add SIR propensities
ReactionNetworkSIR = R6::R6Class(
  classname = "ReactionNetworkSIR",
  inherit = ReactionNetwork,
  public = list(
    get_propensity = function(i,s){
      self$react_rates[i] * switch (i,
        prod(s[1L:2L]), # infection
        s[2L],          # recovery
        1               # arrival of new susceptible
      )
    }
  )
)
# extend to add specific prior
MHSamplerITSAdaptSTGHS2017SIR = R6::R6Class(
  classname = "MHSamplerITSAdaptSTGHS2017SIR",
  inherit = MHSamplerReactNetITSAdaptST,
  public = list(
    # prior logdensity used by the authors
    ldprior = function(theta){sum(dgamma(x=theta,1.5,5,log=TRUE))}
  )
)

# initialize sampler
get_sampler_GHS2017_SIR = function(...){
  updates = matrix(as.integer(c(-1, 1, 0,
                                 0,-1, 1,
                                 1, 0, 0)),nrow=3L,byrow = TRUE)
  theta_true = c(0.4, 0.5, 0.4) # https://github.com/ageorgou/roulette/blob/1625f7dc3b8b55845d1ac8458db9d4418dcdd323/experiments/SIRinfMHServer.m#L10
  par_names = c("Infection", "Recovery", "Arrival")
  CTMC=ReactionNetworkSIR$new( # instantiate a CTMC object
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "obsSIRinf", package = "ctmc3"))
  # dta = load_exp_data(path_to_data = file.path(".","data","obsSIRinf"))
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,theta_true=theta_true)
  return(list(sampler_constructor=MHSamplerITSAdaptSTGHS2017SIR,
              pars=list_pars))
}

###############################################################################
# SG2019 experiments
###############################################################################

#######################################
# Schlögel
#######################################

# extend ReactionNetwork to add Schlögel propensities
ReactionNetworkSch = R6::R6Class(
  classname = "ReactionNetworkSch",
  inherit = ReactionNetwork,
  public = list(
    get_propensity = function(i,s){
      rr=self$react_rates
      switch (i,
        rr[2L]*max(0,s*(s-1)*(s-2)/6) + rr[4L]*s, # Death: 3X -> 2X + A AND X -> B
        rr[1L]*max(0,s*(s-1)/2) + rr[3L]          # Birth: A + 2X -> 3X AND B -> X
      )
    }
  )
)

# prior logdensity used by the authors
# update: changed from normal to lognormal to actually match authors
SG2019_Sch_ldprior = function(theta) sum(dlnorm(x=theta, log=TRUE))

# extend to add specific prior
MHSamplerITSAdaptSTSG2019Sch = R6::R6Class(
  classname = "MHSamplerITSAdaptSTSG2019Sch",
  inherit = MHSamplerReactNetITSAdaptST,
  public = list(ldprior = SG2019_Sch_ldprior))

# extend to add specific prior
MHSamplerRTSAdaptSTSG2019Sch = R6::R6Class(
  classname = "MHSamplerRTSAdaptSTSG2019Sch",
  inherit = MHSamplerReactNetRTSAdaptST,
  public = list(ldprior = SG2019_Sch_ldprior))

# initialize sampler
get_sampler_SG2019_Sch = function(reg_ts=FALSE,...){
  # species: X only, birth-death version of the model (collapsed)
  # like Vellela & Qian (2009), Eqs 2.6,2.7, but with functional forms in Example 2 SG2019
  updates = matrix(c(-1L,1L), nrow=2L, ncol=1L)
  theta_true=c(3, 0.5, 0.5, 3) # Table 5 SG2019
  par_names = c("A + 2X -> 3X","3X -> 2X + A", "B -> X","X -> B")
  CTMC=ReactionNetworkSch$new( # instantiate a CTMC object
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "SG2019_Sch", package = "ctmc3"))
  # dta = load_exp_data(file.path(".","data","SG2019_Sch"))
  sampler_constructor=if(reg_ts) MHSamplerRTSAdaptSTSG2019Sch else MHSamplerITSAdaptSTSG2019Sch
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,
                 theta_true=theta_true)
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}

#######################################
# LV20
# Note: 1 reason why using directions = rows of update matrix doesn't work well here is that in this
# version of LV the only way to expand in both (pred,prey) directions is via
# combinations of 3rd and 2nd reactions. This is much slower than applying the
# 2 independent reactions for "birth" that the LV model above has
# Graphically: you want to reach (100,100) but can only do so by 2 moves:
#     - moving to the right via (0,1)
#     - moving up and left via (1,-1)
# Worst case: need to reach (0,200) and then apply (1,-1) 100 times!!
#######################################

# extend ReactionNetwork to add LV (3 reactions) propensities
ReactionNetworkLV3R = R6::R6Class(
  classname = "ReactionNetworkLV3R",
  inherit = ReactionNetwork,
  public = list(
    get_propensity = function(i,s){
      self$react_rates[i] * switch (i,
        s[1L],   # death of predator
        s[2L],   # birth of prey
        prod(s)  # birth predator + death of prey
      )
    }
  )
)

# prior logdensity used by the authors
# update: changed from normal to lognormal to actually match authors
SG2019_LV_ldprior = function(theta,ml=log(c(0.2,0.2,0.02))){
  sum(dlnorm(x=theta, meanlog = ml, log=TRUE))
}

# extend to add specific prior
MHSamplerITSAdaptSTSG2019LV = R6::R6Class(
  classname = "MHSamplerITSAdaptSTSG2019LV",
  inherit = MHSamplerReactNetITSAdaptST,
  public = list(ldprior = SG2019_LV_ldprior))

# extend to add specific prior
MHSamplerRTSAdaptSTSG2019LV = R6::R6Class(
  classname = "MHSamplerRTSAdaptSTSG2019LV",
  inherit = MHSamplerReactNetRTSAdaptST,
  public = list(ldprior = SG2019_LV_ldprior))

# initialize sampler
get_sampler_SG2019_LV = function(reg_ts=FALSE,...){
  # Lotka-Volterra with 3 reactions [predator:col1 - prey:col2]
  updates = matrix(as.integer(c(-1, 0,  # death of predator
                                 0, 1,  # birth of prey
                                 1,-1)),# birth pred + death prey
                   nrow=3L,byrow = TRUE)
  theta_true=c(0.3, 0.4, 0.01) # table 5 in SG2019
  par_names = c("Death of predator","Birth of prey","Pred+Prey -> 2Pred")
  CTMC=ReactionNetworkLV3R$new( # instantiate a CTMC object
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "SG2019_LV20", package = "ctmc3"))
  # dta = load_exp_data(file.path(".","data","SG2019_LV20"))
  sampler_constructor=if(reg_ts) MHSamplerRTSAdaptSTSG2019LV else MHSamplerITSAdaptSTSG2019LV
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,
                 theta_true=theta_true)
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}

##############################################################################
# other models not used in experiments
##############################################################################

#######################################
# M/M/c queue => uniformizable
#######################################

# extend ReactionNetwork to add M/M/c propensities
ReactionNetworkMMc = R6::R6Class(
  classname = "ReactionNetworkMMc",
  inherit = ReactionNetwork,
  public = list(
    nserv = NULL, # number of servers (c)
    initialize = function(nserv,...) {
      self$nserv=nserv
      super$initialize(...)
    },
    get_propensity = function(i,s){
      switch (i,
              self$react_rates[1L]*min(s,self$nserv), # Death
              self$react_rates[2L]                    # Birth
      )
    }
  )
)

# initialize sampler
get_sampler_MMc = function(reg_ts=FALSE,nserv=5L,...){
  # 1 species: number of customers being served
  updates = matrix(c(-1L,1L), nrow=2L, ncol=1L)
  theta_true=c(2, 8) # assume nserv is known
  par_names = c("Service rate", "Arrival rate")
  CTMC=ReactionNetworkMMc$new( # instantiate a CTMC object
    nserv        = nserv,
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "MMc", package = "ctmc3"))
  # dta = load_exp_data(file.path(".","data","MMc"))
  sampler_constructor=if(reg_ts) MHSamplerRTSAdaptSTSG2019Sch else MHSamplerITSAdaptSTSG2019Sch
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,
                 theta_true=theta_true)
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}

