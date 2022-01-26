###############################################################################
# define classes
###############################################################################

#######################################
# Infinite CTMC (basic required methods)
#######################################

InfiniteCTMC = R6::R6Class(
  classname = "InfiniteCTMC",
  public = list(
    simulate = function(...){stop("not implemented")}, # simulate (fully observed) process
    sim_data = function(s0,Tend=1,n_obs=25L,min_jumps=0L,
                        sim_trace,sample_times){
      if(missing(sim_trace)){
        sim_trace = self$simulate(s0=s0,Tend=Tend) # simulate process
        while(length(sim_trace$t) <= min_jumps){ # reject if too few jumps
          sim_trace = self$simulate(s0=s0,Tend=Tend)
        }
      }

      if(missing(sample_times)){
        # stratified sampling
        # sample 3 times inside each interval to allow for repeatedly sampling the same state
        sample_times = numeric(0L)
        for(i in seq.int(2L,length(sim_trace$t))){
          sample_times = c(
            sample_times,
            runif(n = 3L,#n = ceiling(max(sim_trace$t)/(sim_trace$t[i]-sim_trace$t[i-1L])),
                  min = sim_trace$t[i-1L],max = sim_trace$t[i])
          )
        }
        sample_times = sort(sample(sample_times, n_obs))
      }

      # build subsample object
      dta = list(
        t = c(0, sample_times),
        # for each sampled time, find closest jump time to the left, and assign that state
        s = unname(rbind(s0, sim_trace$s[
          vapply(sample_times, function(u) {
            d=u - sim_trace$t
            which.min(ifelse(d>=0,d,Inf))
          }, FUN.VALUE = integer(1L)),,drop=FALSE]))
      )
      dta$t = matrix(dta$t, ncol = 1L)
      dta$dt = rbind(NA_real_,diff(dta$t)) # compute delta t
      dta$s_pre = rbind(rep(NA_integer_,ncol(sim_trace$s)),
                        dta$s[-nrow(dta$s),,drop=FALSE]) # add lagged state
      # dta = lapply(dta, function(l){l[-1L,,drop=FALSE]}) # can now drop observation at t=0
      return(list(dta=dta,sim_trace=sim_trace))
    }
  )
)

##############################################################################
# Reaction Networks
##############################################################################

ReactionNetwork = R6::R6Class(
  classname = "ReactionNetwork",
  inherit = InfiniteCTMC,
  public = list(
    updates = NULL, # (int matrix) nrow=n_react,ncol=n_spec
    react_rates = NULL, # (>0 vector) length n_react
    initialize = function(updates,react_rates) {
      self$updates=updates
      self$react_rates=react_rates
    },

    # evaluate propensity of a rate at a given state (must be implemented by extending class)
    get_propensity = function(...){stop("not implemented")},
    get_all_propensities=function(s){
      vapply(seq_len(nrow(self$updates)), function(i){
        self$get_propensity(i,s)
      }, FUN.VALUE = numeric(1L))
    },

    # simulate a sample path
    simulate = function(s0, Tend, max_jumps=1000L){
      t = 0 # set initial time
      t_vec = t # initialize jump times vector
      s_mat = matrix(s0,ncol=ncol(self$updates)) # initialize state
      i = 0L # row tracker in s_mat
      while(i < max_jumps && t<Tend){
        i = i+1L
        clock_rates = self$get_all_propensities(s_mat[i,,drop=FALSE])
        react_times = suppressWarnings(rexp(n=length(clock_rates), rate = clock_rates))
        react_times = ifelse(is.nan(react_times), Inf, react_times)
        r = which.min(react_times)
        s_mat = rbind(s_mat, s_mat[i,,drop=FALSE] + self$updates[r,,drop=FALSE])
        t = t + react_times[r]
        t_vec = c(t_vec, t)
      }
      if(t>Tend){
        s_mat = s_mat[-(i+1L),,drop=FALSE]
        t_vec = t_vec[-(i+1L)]
      }
      return(list(s=s_mat, t=t_vec))
    }
  )
)
