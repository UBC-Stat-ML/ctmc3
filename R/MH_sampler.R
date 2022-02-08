###############################################################################
# Generic Metropolis-Hastings sampler
###############################################################################

MHSampler = R6::R6Class(
  classname = "MHSampler",
  public = list(
    theta_0    = NULL, # initial value of parameter for the chain
    dim        = NULL, # dimension of space
    theta      = NULL, # matrix S x dim with values of the sampled chain
    theta_true = NULL, # true value of theta if available
    ld_stats   = NULL, # matrix with trace of logdensities (prior, loglik, target)
    log_file   = NULL, # file connection for logging
    par_names  = NULL, # optional names for theta
    initialize = function(theta_0,log_file=NULL,par_names=NULL,theta_true=NULL) {
      self$theta_0    = theta_0
      self$dim        = length(self$theta_0)
      self$log_file   = if(is.null(log_file)) nullfile() else log_file
      self$par_names  = par_names
      self$theta_true = theta_true
    },
    rprop     = function(...){stop("not implemented")}, # RNG for proposal kernel
    ldprop    = function(...){stop("not implemented")}, # log-density of proposal kernel
    ldtarget  = function(...){stop("not implemented")}, # returns list with log-density of "target" distribution (mandatory), "loglik" and "prior" (both optional; i.e., possibly NA)
    run_chain = function(S=1000L,print_every=100L){
      # allocate storage and initialize
      theta              = matrix(NA_real_, nrow = S, ncol = self$dim)
      colnames(theta)    = names(self$theta_0)
      theta[1L,]         = self$theta_0                             # init chain
      res                = self$ldtarget(self$theta_0)              # evaluate target at starting point
      ld_stats           = matrix(NA_real_,nrow=S,ncol=length(res)) # define storage for trace of ldtarget
      colnames(ld_stats) = names(res)                               # set names
      ld_stats[1L,]      = unlist(res)                              # store current
      cur_ld             = res$target                               # cache current logdensity of target
      
      if(is.infinite(cur_ld) && cur_ld<0) stop("log-target(theta_0)=-Inf, cannot proceed.")
      
      # sampling loop
      n_reject = 0L # counter of rejections
      for(s in seq.int(2L,S)){
        if(print_every>0L && (s %% print_every == 0L))
          cat(sprintf("Step %d, rejection rate so far %.0f%%\n", s, round(100*n_reject/(s-1))))
        ptheta  = self$rprop(theta[s-1L,]) # draw proposal
        res     = self$ldtarget(ptheta)    # compute target at proposal (e.g. for Bayesian model involves computing prior and loglik)
        prop_ld = res$target               # get target logdensity at proposal
        laccept = (prop_ld + self$ldprop(cur=ptheta, prop=theta[s-1L,])) -
          (cur_ld + self$ldprop(cur=theta[s-1L,],prop=ptheta))
        if(rexp(1L) + laccept >= 0){
          theta[s,]    = ptheta  # move!
          cur_ld       = prop_ld # update cache of logdensity
          ld_stats[s,] = unlist(res)
        } else{
          n_reject     = n_reject+1L
          theta[s,]    = theta[s-1L,]
          ld_stats[s,] = ld_stats[s-1L,]
        }
      }
      self$theta    = theta
      self$ld_stats = ld_stats
      stats         = list(n_reject=n_reject)
      invisible(
        list(
          theta = `colnames<-`(self$theta,self$par_names),
          ld_stats=self$ld_stats, stats = stats)
      )
    },

    # run in parallel using forking
    run_chains_parallel = function(S,n_chains,print_every=100L){
      fproc=lapply(seq_len(n_chains), function(j){
        parallel::mcparallel(expr = {
          sink(file = self$log_file, append = TRUE)
          res = self$run_chain(S,print_every=print_every)
          sink()
          res
        })
      })
      res=parallel::mccollect(fproc)
      ind_err=which(sapply(res, is, class2 = "try-error"))
      if(length(ind_err)>0L) stop(res[[ind_err[1L]]]) # break on first error
      return(res)
    },

    # set theta_0 close to a mode estimated using output from run_chains_parallel
    set_theta_0_from_par_run = function(res_run){
      # find max within chains
      ind_vec = sapply(res_run, function(l) {
        which.max(l$ld_stats[, "target"])
      })
      # find max across chains
      ind = which.max(sapply(seq_along(res_run), function(i) {
        res_run[[i]]$ld_stats[ind_vec[i], "target"]
      }))
      self$theta_0 = res_run[[ind]]$theta[ind_vec[ind], ]
      list(theta=self$theta_0,
           target=res_run[[ind]]$ld_stats[ind_vec[ind], "target"])
    }
  )
)
