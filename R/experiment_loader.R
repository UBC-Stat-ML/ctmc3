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
get_sampler = function(
  exp_name,
  reg_ts,
  gtp_solver,
  tuned_pars_path,
  tuned_pars      = TRUE,
  use_theta_true  = TRUE,
  ...)
  {
  caller = switch(
    EXPR = exp_name,
    GHS2017_LV  = get_sampler_GHS2017_LV,
    GHS2017_SIR = get_sampler_GHS2017_SIR,
    SG2019_LV20 = get_sampler_SG2019_LV,
    SG2019_Sch  = get_sampler_SG2019_Sch,
    MMc         = get_sampler_MMc,
    SG2019_Sch_log  = get_sampler_SG2019_Sch_log,
    SG2019_LV20_log = get_sampler_SG2019_LV_log,
    GHS2017_LV_log  = get_sampler_GHS2017_LV_log,
    stop("Unknown experiment."))
  res_list = do.call(caller,list(reg_ts=reg_ts)) # get sampler constructor + experiment parameters
  
  # compile list of arguments to pass to the constructor
  alist=c(list(exp_name=exp_name,gtp_solver=gtp_solver),res_list$pars,list(...))
  if(!missing(tuned_pars_path) || tuned_pars){                   # use stored parameters
    model_string = paste0(
      sprintf("%s_%s%s",exp_name,ifelse(reg_ts,"reg_ts_",""), gtp_solver)
    )
    # build paths to files
    if(missing(tuned_pars_path))
      tuned_pars_path = file.path(
        system.file(package="ctmc3"), "extdata", "MH_sampler_tuning"
      )
    model_pars     = file.path(tuned_pars_path, model_string)
    path_to_varmat = file.path(model_pars, "varmat.txt")
    path_to_theta0 = file.path(model_pars, "theta_0.txt")
    path_to_st_dta = file.path(model_pars, "st_dta.txt")
    
    # load parameters
    alist$varmat = load_varmat(
      path_to_varmat = path_to_varmat,
      par_names      = alist$par_names
    )
    alist$theta_0        = scan(path_to_theta0, quiet = TRUE)
    names(alist$theta_0) = alist$par_names
    st_dta               = read.table(path_to_st_dta, header = TRUE)
    
    # build sampler and set stopping times
    sampler = do.call("new", alist, envir = res_list$sampler_constructor)
    sampler$st_list = lapply(
      seq_len(nrow(st_dta)),
      function(i){ do.call(GeometricST$new, st_dta[i,]) } # note we use here a dataframe row as list
    )
  }else{
    alist$theta_0 = if(use_theta_true) alist$theta_true else get_rnd_theta(alist$theta_true)
    alist$varmat  = diag(ifelse(grepl("_log$", exp_name),0.01,0.1)*alist$theta_0^2)
    sampler       = do.call("new", alist, envir = res_list$sampler_constructor)
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

# prior logdensity used by the authors
GHS2017_LV_ldprior = function(theta){
  sum(dgamma(x = theta, 4, 10000, log = TRUE)) # == -Inf if any(theta<0)
}
# normal approximation to the log scale prior (used for comparison to SG19)
# use simple moment matching, by leveraging that for G ~ Gamma(a,b), then
# - E[log(G)]   = digamma(a) - log(b)
# - Var[log(G)] = trigamma(a)
GHS2017_LV_log_ldprior = function(theta,m=digamma(4)-log(10000),s=sqrt(trigamma(4))){
  sum(dnorm(x = theta, mean = m, sd = s, log = TRUE))
}

# extend to add specific prior
MHSamplerITSGHS2017LV = R6::R6Class(
  classname = "MHSamplerITSGHS2017LV",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = GHS2017_LV_ldprior))
MHSamplerRTSGHS2017LV = R6::R6Class(
  classname = "MHSamplerRTSGHS2017LV",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = GHS2017_LV_ldprior))
MHSamplerITSGHS2017LVLog = R6::R6Class(
  classname = "MHSamplerITSGHS2017LVLog",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = GHS2017_LV_log_ldprior))
MHSamplerRTSGHS2017LVLog = R6::R6Class(
  classname = "MHSamplerRTSGHS2017LVLog",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = GHS2017_LV_log_ldprior))

# initialize sampler
get_sampler_GHS2017_LV = function(reg_ts=FALSE,use_log=FALSE,...){
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
  sampler_constructor = if(reg_ts){
    if(use_log) MHSamplerRTSGHS2017LVLog else MHSamplerRTSGHS2017LV 
  }else{
    if(use_log) MHSamplerITSGHS2017LVLog else MHSamplerITSGHS2017LV
  }
  list_pars=list(
    CTMC=CTMC,dta=dta,par_names=par_names,
    theta_true=if(use_log) log(theta_true) else theta_true,
    theta_to_rrs=if(use_log) function(x){exp(x)} else function(x){x}
  )
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}
get_sampler_GHS2017_LV_log = function(...) get_sampler_GHS2017_LV(use_log = TRUE, ...)

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
MHSamplerITSGHS2017SIR = R6::R6Class(
  classname = "MHSamplerITSGHS2017SIR",
  inherit = MHSamplerReactNetITS,
  public = list(
    # prior logdensity used by the authors
    ldprior = function(theta){sum(dgamma(x=theta,1.5,5,log=TRUE))} # == -Inf if any(theta<0)
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
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,theta_true=theta_true)
  return(list(sampler_constructor=MHSamplerITSGHS2017SIR,
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
SG2019_Sch_ldprior = function(theta) sum(dlnorm(x=theta, log=TRUE))    # == -Inf if any(theta<0)
SG2019_Sch_log_ldprior = function(theta) sum(dnorm(x=theta, log=TRUE)) # work in log space

# extend classes to add specific prior
MHSamplerITSSG2019Sch = R6::R6Class(
  classname = "MHSamplerITSSG2019Sch",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = SG2019_Sch_ldprior))
MHSamplerRTSSG2019Sch = R6::R6Class(
  classname = "MHSamplerRTSSG2019Sch",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = SG2019_Sch_ldprior))
MHSamplerITSSG2019SchLog = R6::R6Class(
  classname = "MHSamplerITSSG2019SchLog",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = SG2019_Sch_log_ldprior))
MHSamplerRTSSG2019SchLog = R6::R6Class(
  classname = "MHSamplerRTSSG2019SchLog",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = SG2019_Sch_log_ldprior))

# initialize sampler
get_sampler_SG2019_Sch = function(reg_ts = FALSE, use_log = FALSE, ...){
  # species: X only, birth-death version of the model (collapsed)
  # like Vellela & Qian (2009), Eqs 2.6,2.7, but with functional forms in Example 2 SG2019
  updates = matrix(c(-1L,1L), nrow=2L, ncol=1L)
  theta_true=c(3, 0.5, 0.5, 3) # Table 5 SG2019
  par_names = c("A + 2X -> 3X","3X -> 2X + A", "B -> X","X -> B")
  CTMC=ReactionNetworkSch$new( # instantiate a CTMC object
    updates      = updates,
    react_rates  = theta_true)
  dta = load_exp_data(system.file("extdata", "SG2019_Sch", package = "ctmc3"))
  sampler_constructor=if(reg_ts){
    if(use_log) MHSamplerRTSSG2019SchLog else MHSamplerRTSSG2019Sch
  }else{
    if(use_log) MHSamplerITSSG2019SchLog else MHSamplerITSSG2019Sch
  }
  list_pars=list(
    CTMC=CTMC,dta=dta,par_names=par_names,
    theta_true=if(use_log) log(theta_true) else theta_true,
    theta_to_rrs=if(use_log) function(x){exp(x)} else function(x){x}
  )
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}
get_sampler_SG2019_Sch_log = function(...) get_sampler_SG2019_Sch(use_log = TRUE, ...)

#######################################
# LV20
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
SG2019_LV_ldprior = function(theta, ml = log(c(0.2,0.2,0.02))){
  sum(dlnorm(x=theta, meanlog = ml, log=TRUE))
}
SG2019_LV_log_ldprior = function(theta, m = log(c(0.2,0.2,0.02))){
  sum(dnorm(x=theta, mean = m, log=TRUE))
}

# extend to add specific prior
MHSamplerITSSG2019LV = R6::R6Class(
  classname = "MHSamplerITSSG2019LV",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = SG2019_LV_ldprior))
MHSamplerRTSSG2019LV = R6::R6Class(
  classname = "MHSamplerRTSSG2019LV",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = SG2019_LV_ldprior))
MHSamplerITSSG2019LVLog = R6::R6Class(
  classname = "MHSamplerITSSG2019LVLog",
  inherit = MHSamplerReactNetITS,
  public = list(ldprior = SG2019_LV_log_ldprior))
MHSamplerRTSSG2019LVLog = R6::R6Class(
  classname = "MHSamplerRTSSG2019LVLog",
  inherit = MHSamplerReactNetRTS,
  public = list(ldprior = SG2019_LV_log_ldprior))

# initialize sampler
get_sampler_SG2019_LV = function(reg_ts = FALSE, use_log = FALSE, ...){
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
  sampler_constructor=if(reg_ts){
    if(use_log) MHSamplerRTSSG2019LVLog else MHSamplerRTSSG2019LV 
  }else{
    if(use_log) MHSamplerITSSG2019LVLog else MHSamplerITSSG2019LV
  }
  list_pars=list(
    CTMC=CTMC,dta=dta,par_names=par_names,
    theta_true=if(use_log) log(theta_true) else theta_true,
    theta_to_rrs=if(use_log) function(x){exp(x)} else function(x){x}
  )
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}
get_sampler_SG2019_LV_log = function(...) get_sampler_SG2019_LV(use_log = TRUE, ...)
                                                                
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
  sampler_constructor=if(reg_ts) MHSamplerRTSSG2019Sch else MHSamplerITSSG2019Sch
  list_pars=list(CTMC=CTMC,dta=dta,par_names=par_names,
                 theta_true=theta_true)
  return(list(sampler_constructor=sampler_constructor,pars=list_pars))
}

