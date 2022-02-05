###############################################################################
# Extend Irregular TS sampler to add adaptive stop times and unbiased estimator
###############################################################################

MHSamplerReactNetITSAdaptST = R6::R6Class(
  classname = "MHSamplerReactNetITSAdaptST",
  inherit = MHSamplerReactNetITS,
  public = list(
    st_list=NULL, # list of adaptive stop-times, one per observation
    dta_adapt=NULL, # info about convergence of approx, one per observation

    ######################################################################
    # methods for stopping times
    ######################################################################

    study_convergence = function(theta=self$theta_0,
                                 max_err_thresh = sqrt(.Machine$double.eps),
                                 N_lo_eps = if(self$gtp_solver=="skeletoid") -10L else 0L){
      self$set_params(theta) # translate theta to CTMC's parameters, and set it
      dta_adapt = vector("list", self$n_obs)
      for(ind_obs in seq_len(self$n_obs)){
        # ind_obs=10L
        cat(sprintf("\nStudying convergence for obs %d\n",ind_obs))

        # set N_eps high and study truncation Cauchy error
        # note: first Cauchy error is undefined so set arbitrary
        # note: need to distinguish between low error and 0 probability of jump
        # if the data is correct, then the latter occurs only because approximation
        # is too poor to admit jumps that involve many intermediate steps
        # Because of the above, a "bump" appears in the profile of the Cauchy error
        # We try to detect being in left tail with heuristics (see below)
        N_hi_eps=as.integer(-log10(max_err_thresh))
        N_trunc=0L; old_est=0; cauchy_err=1; dta_trunc=data.frame()
        prob_vec=numeric() # storage for estimates
        exploring_left_tail = TRUE # still on the left tail of the bump
        while(exploring_left_tail || cauchy_err>max_err_thresh){
          new_est = self$trans_prob_est(
            o=ind_obs,N_trunc = N_trunc,N_eps = N_hi_eps)
          prob_vec=c(prob_vec,new_est)
          cauchy_err = new_est-old_est
          old_est=new_est
          dta_trunc=rbind(
            dta_trunc, data.frame(N=N_trunc,cauchy_err=cauchy_err))
          if(self$debug)cat(sprintf("N_trunc=%d, prob=%g, Cauchy-err=%g\n",
                                    N_trunc,new_est,cauchy_err))
          N_trunc=N_trunc+1L
          exploring_left_tail=
            (sum(dta_trunc$cauchy_err<=100*.Machine$double.eps)<=8L && # This controls how long we explore before giving up: <=0 Cauchy err are numerical errors that occur when solver is as exact as possible. If we see "too many" it means there's nothing else to do
               (diff(range(prob_vec)) < max_err_thresh || # estimated trans-prob is fixed at the first estimate with no improvement
                  all(prob_vec<max_err_thresh))) # all estimated trans-prob are tiny
        }
        N_hi_trunc = N_trunc-1L # this truncation level gives Cauchy err < max_err_thresh
        best_guess = max(prob_vec)
        dta_trunc$cauchy_err[1L]=NA_real_
        dta_trunc$err = best_guess-prob_vec
        dta_trunc$prob_vec = prob_vec
        # delete unsolicited exploration of right tail (because we detected too late we passed the last peak)
        ind_enough=max(2L,which(abs(prob_vec-best_guess)<max_err_thresh)[1L])
        dta_trunc=dta_trunc[seq_len(ind_enough),]

        # set N_eps high and study solver approximation error (have best guess)
        N_eps=N_lo_eps; err=1; dta_eps=data.frame(); prob_vec=numeric()
        while(abs(err)>max_err_thresh){
          new_est = self$trans_prob_est(
            o=ind_obs,N_trunc = N_hi_trunc,N_eps = N_eps)
          prob_vec = c(prob_vec,new_est)
          err = best_guess-new_est
          dta_eps=rbind(dta_eps, data.frame(N=N_eps,err=err))
          N_eps=N_eps+1L
        }
        dta_eps$prob_vec=prob_vec
        dta_eps$cauchy_err=c(NA_real_,diff(dta_eps$prob_vec))
        ind_rep = which(dta_eps$err==dta_eps$err[1L]) # remove duplicates in the left end
        if(length(ind_rep)>1L) dta_eps=dta_eps[-seq.int(1L,length(ind_rep)-1L),]

        dta_adapt[[ind_obs]]=list(dta_trunc=dta_trunc,dta_eps=dta_eps)
      }
      self$dta_adapt=dta_adapt
    },

    # to set the offset, we base our decision on the optimal stopping pmf
    # for the estimator
    #       Z=D_n/p(n)
    # where D_0=X_O and D_n=X_n-X_{n-1}. The optimal 0 variance pmf is
    #       pmf(n) = D_n/X \propto D_n
    # with X := sum_{n>=0} D_n.
    # (note: this is not the the estimator we end up using, but
    #  this is still a good heuristic for setting O)
    # Therefore, we set offset := min{n: X_n/X >= min_mass}
    # However, we also need pmf(n) decreasing after the offset. Recall that
    # the estimator we use is
    #       Z = X_O + (X_N+1 - X_N)/pmf(N)
    # and its optimal pmf is
    #       pmf(n) = (X_n+1 - X_n)/(X - X_O) \propto (X_n+1 - X_n)
    # Therefore, if we want to mimic this with a simple Geom stop time,
    # we must ensure (X_n+1 - X_n) is decreasing
    find_offset = function(dta_conv,min_mass,
                           max_mass_post_offset=100*.Machine$double.eps,
                           meaningful_decrease_cerr=1E-6){
      # dta_conv=self$dta_adapt[[2L]]$dta_trunc
      # compute pmf0, and filter data to discard values too far to the right
      pmf_0=dta_conv$prob_vec/max(dta_conv$prob_vec)
      ind_last=which(pmf_0>1-max_mass_post_offset)[1L]
      tryCatch(
        {dta_conv=dta_conv[seq_len(ind_last),]}
        , error=function(e){
          cat("\nERROR! Dumping pmf_0, then dta_conv\n\n")
          print(pmf_0);print(dta_conv);stop(e)
        })
      if(nrow(dta_conv)==1L) return(1L)

      # find last peak to ensure monotone decreasing
      if(nrow(dta_conv)==2L){
        ind_origin=1L # not enough info to check for peaks (first Cauchy err is always NA)
      }else{ # have at least 2 non-NA Cauchy err's
        shift_cerr=dta_conv$cauchy_err[-1L]
        # plot(shift_cerr,log="y");lines(shift_cerr)
        drev=diff(rev(shift_cerr))
        # plot(drev);lines(drev)
        ind_flip_sign=which(drev< -meaningful_decrease_cerr)[1L]
        ind_origin=length(shift_cerr)-ind_flip_sign+1L
        if(is.na(ind_origin)) ind_origin=1L # already monotone decreasing
      }

      # find offset via pmf0 criterion, correct with last peak info
      ind_origin=max(ind_origin,which(pmf_0>min_mass)[1L])
      return(ind_origin) # note: returns index not an actual N
    },
    
    # called only within set_AST (below)
    study_joint_convergence_one_obs = function(
      ind_obs,
      min_mass      = 0.8,
      min_eps_speed = if(self$correct_unif) 20.0 else 0.1,
      tol           = sqrt(.Machine$double.eps)
    ){
      # find preliminary offsets and eps_speed, and init storage
      dta_adapt     = self$dta_adapt[[ind_obs]]
      O_trunc       = with(dta_adapt, dta_trunc$N[self$find_offset(dta_trunc, min_mass)])
      O_eps         = with(dta_adapt, dta_eps$N[self$find_offset(dta_eps, min_mass)])
      N_hi_trunc    = max(dta_adapt$dta_trunc$N)
      N_hi_eps      = max(dta_adapt$dta_eps$N)
      N_trunc_vec   = seq.int(O_trunc,N_hi_trunc)
      prob_vec      = numeric(length(N_trunc_vec))
      eps_speed     = max(min_eps_speed,(N_hi_eps-O_eps)/length(N_trunc_vec)) # init at default
      N_eps_vec     = numeric(length(N_trunc_vec))
      N_eps_vec[1L] = O_eps
      prob_vec[1L]  = self$trans_prob_est(ind_obs,N_trunc = O_trunc,N_eps = O_eps)
      if(self$debug) cat(sprintf(", O_trunc_pre=%d, O_eps_pre=%.1f",O_trunc,O_eps))
      
      # explore joint convergence, check monotonicity holds (it must except for unif)
      i             = 1L
      N_trunc       = N_trunc_vec[i]
      N_eps         = O_eps
      ind_last_doub = 1L
      prob          = prob_vec[1L]
      prob_true     = max(dta_adapt$dta_trunc$prob_vec, dta_adapt$dta_eps$prob_vec)
      while(N_trunc<N_hi_trunc && prob+tol<prob_true){
        i       = i+1L
        N_trunc = N_trunc_vec[i]
        N_eps   = N_eps_vec[i-1L]+eps_speed
        prob    = self$trans_prob_est(ind_obs, N_trunc = N_trunc, N_eps = N_eps)
        
        # this loop corrects for the potential non-double-monotonicity of unif
        # by increasing eps_speed until the joint sequence is increasing
        while(prob-prob_vec[i-1L]<0 && prob+tol<prob_true){
          if(self$debug) cat(", eps_speed doubled")
          eps_speed     = 2*eps_speed
          ind_last_doub = i
          N_eps         = N_eps_vec[i-1L]+eps_speed
          prob          = self$trans_prob_est(ind_obs,N_trunc = N_trunc,N_eps = N_eps)
        }
        
        N_eps_vec[i] = N_eps
        prob_vec[i]  = prob
      }
      if(self$debug) cat(sprintf(", eps_speed=%.1f", eps_speed))
      
      # build the joint convergence dataframe
      # if we doubled after the first loop step, we need to reconstruct prov_vec
      N_trunc_vec = N_trunc_vec[1L:i]
      if(ind_last_doub <= 2L){
        N_eps_vec = N_eps_vec[1L:i]
        prob_vec  = prob_vec[1L:i]
      }else{                                       # need to reconstruct joint seq
        N_eps_vec        = seq(O_eps,by=eps_speed,length.out = i)
        prob_vec_new     = numeric(i)
        prob_vec_new[1L] = prob_vec[1L]              # at least the origin is correct
        for(j in seq.int(2L,i)){
          prob_vec_new[j] = self$trans_prob_est(
            ind_obs, N_trunc = N_trunc_vec[j], N_eps = N_eps_vec[j]
          )
        }
        prob_vec = prob_vec_new
      }
      dta_joint=data.frame(                        # store results
        N          = N_trunc_vec,
        N_trunc    = N_trunc_vec, 
        N_eps      = N_eps_vec,
        prob_vec   = prob_vec, 
        cauchy_err = c(NA_real_,diff(prob_vec))
      )
      
      # store results and return eps_speed
      self$dta_adapt[[ind_obs]]$dta_joint = dta_joint
      return(eps_speed)
    },
    
    set_AST = function(
      min_mass    = 0.8,
      slope_alpha = 0.99,
      min_p_geom  = 0.1,
      max_p_geom  = 0.9
      ){
      for(ind_obs in seq_len(self$n_obs)){
        # ind_obs=5L
        if(self$debug) cat(sprintf("set_AST: obs=%d",ind_obs))
        
        # study joint convergence and set eps_speed
        # note: joint convergence data is built using a given min_mass so it cannot
        # be reused. it must be recalculated every time we change min_mass
        eps_speed = self$study_joint_convergence_one_obs(
          ind_obs = ind_obs, min_mass = min_mass
        )
    
        # find joint offsets
        dta_joint  = self$dta_adapt[[ind_obs]]$dta_joint
        ind_origin = self$find_offset(dta_joint, min_mass)
        O_trunc    = dta_joint$N_trunc[ind_origin]
        O_eps      = dta_joint$N_eps[ind_origin]
        if(self$debug) cat(sprintf(", O_trunc=%d, O_eps=%.1f",O_trunc,O_eps))
        
        # set parameter of geometric stop time by mimicking optimal stopping
        # probability: pmf(n) \propto X_n+1 - X_n
        if(ind_origin>=nrow(dta_joint)-3L){          # check if we have enough data for regression (>=5 points)
          p_geom=max_p_geom                          # we're at the very end of acceptable error
        }else{
          shft=min(0,min(dta_joint$cauchy_err[-1L]))
          stop_prob=dta_joint$cauchy_err[-1L]-shft # shift to ensure >=0
          stop_prob=if(length(stop_prob)==1L) 1 else stop_prob/sum(stop_prob)
          dta_joint$stop_prob=c(stop_prob,0)
          # with(dta_joint[-nrow(dta_joint),],plot(N,stop_prob,log="y",type = "l"))
          # abline(v=O_trunc,lty="dotted")
          # abline(h=stop_prob[ind_origin],lty="dotted")
          # define subset of dta_joint to use for regression
          # then center variables
          ind_reg=seq.int(ind_origin,length(stop_prob)) # from offset to last point
          ind_reg=ind_reg[which(dta_joint$stop_prob[ind_reg]>0)] # need to remove potential 0 introduced by shifting
          dta_lm=dta_joint[ind_reg,];dta_lm$lsp=log(dta_lm$stop_prob)
          dta_lm$lsp=dta_lm$lsp-dta_lm$lsp[1L]; dta_lm$N=dta_lm$N-dta_lm$N[1L]
          fit_lm=lm(lsp~0+N,data=dta_lm) # force 0 intercept to only capture slope
          if(length(ind_reg)<=20L){
            min_slope=-fit_lm$coefficients["N"]
          }else{
            min_slope=-confint(fit_lm,level = slope_alpha)["N",2L]
          }
          # curve({stop_prob[ind_origin]*exp(-min_slope*(x-O_trunc))},
          #       from=O_trunc,to=max(dta_joint$N_trunc)-1L,add=TRUE,lty="dashed")
          p_geom=-unname(expm1(-min_slope))#=1-exp(-min_slope)
          p_geom=max(min_p_geom,min(max_p_geom,p_geom)) # enforce safety limits
        }
        if(self$debug)cat(sprintf(", p_geom=%.1f. done!\n\n",p_geom))

        # create and store a GeometricST object
        stop_time = GeometricST$new(
          O_trunc   = O_trunc,
          O_eps     = O_eps,
          prob      = p_geom,
          eps_speed = eps_speed
        )
        self$st_list[[ind_obs]] = stop_time
      } # end ind_obs loop
    }
  )
)

###############################################################################
# fin
###############################################################################

###############################################################################
# Extend Irregular TS sampler to add adaptive stop times
###############################################################################

MHSamplerReactNetRTSAdaptST = R6::R6Class(
  classname = "MHSamplerReactNetRTSAdaptST",
  inherit = MHSamplerReactNetRTS,
  public = list(
    st_list=NULL, # list of adaptive stop-times, one per observation
    dta_adapt=NULL, # info about convergence of approx, one per observation

    ######################################################################
    # methods for stopping times
    ######################################################################

    study_convergence = function(theta=self$theta_0,
                                 max_err_thresh = sqrt(.Machine$double.eps),
                                 min_ll_thresh = -700,
                                 N_lo_eps = if(self$gtp_solver=="skeletoid") -10L else 0L){
      cat("\nStudying convergence of the likelihood\n")

      # set N_eps high and study truncation Cauchy error
      N_hi_eps=as.integer(-log10(max_err_thresh))
      N_trunc=0L; old_ll=min_ll_thresh; ll_cauchy_err=1; dta_trunc=data.frame()
      ll_vec=numeric() # storage for estimates
      exploring_left_tail = TRUE # still on the left tail of the bump
      while(exploring_left_tail || ll_cauchy_err>max_err_thresh){
        new_ll = self$loglik_biased(theta,N_trunc,N_hi_eps)
        ll_vec=c(ll_vec,new_ll)
        ll_cauchy_err = new_ll-old_ll
        old_ll=new_ll
        dta_trunc=rbind(
          dta_trunc, data.frame(N=N_trunc,ll_cauchy_err=ll_cauchy_err))
        if(self$debug)cat(sprintf("N_trunc=%d, ll=%g, ll_Cauchy-err=%g\n",
                                  N_trunc, new_ll, ll_cauchy_err))
        N_trunc=N_trunc+1L
        exploring_left_tail=
          (sum(dta_trunc$ll_cauchy_err<=100*.Machine$double.eps)<=8L &&
             (diff(range(ll_vec)) < max_err_thresh ||
                all(ll_vec<min_ll_thresh)))
      }
      N_hi_trunc = N_trunc-1L # this truncation level gives Cauchy err < max_err_thresh
      best_guess = ll_vec[N_hi_trunc+1L]
      dta_trunc$ll_cauchy_err[1L]=NA_real_
      dta_trunc$ll_err = best_guess-ll_vec
      dta_trunc$ll_vec = ll_vec

      # set N_eps high and study solver approximation error (have best guess)
      N_eps=N_lo_eps; ll_err=1; dta_eps=data.frame(); ll_vec=numeric()
      while(ll_err>max_err_thresh){
        new_ll = self$loglik_biased(theta,N_hi_trunc,N_eps)
        ll_vec = c(ll_vec,new_ll)
        ll_err = best_guess-new_ll
        dta_eps=rbind(dta_eps, data.frame(N=N_eps,ll_err=ll_err))
        N_eps=N_eps+1L
      }
      dta_eps$ll_vec=ll_vec
      ind_rep = which(dta_eps$ll_err==dta_eps$ll_err[1L]) # remove duplicates in the left end
      if(length(ind_rep)>1L) dta_eps=dta_eps[-seq.int(1L,length(ind_rep)-1L),]
      dta_eps$ll_cauchy_err=c(NA_real_,diff(dta_eps$ll_vec))
      # with(dta_eps,plot(x=N,y=cauchy_err,type="l",log="y"))
      self$dta_adapt=list(list(dta_trunc=dta_trunc,dta_eps=dta_eps))
    },

    # to set the offset, we base our decision on the optimal stopping pmf
    # for the estimator
    #       Z=D_n/p(n)
    # where D_0=X_O and D_n=X_n-X_{n-1}. The optimal 0 variance pmf is
    #       pmf(n) = D_n/X \propto D_n
    # with X := sum_{n>=0} D_n.
    # (note: this is not the the estimator we end up using, but
    #  this is still a good heuristic for setting O)
    # Therefore, we set offset := min{n: X_n/X >= min_mass}
    # However, we also need pmf(n) decreasing after the offset. Recall that
    # the estimator we use is
    #       Z = X_O + (X_N+1 - X_N)/pmf(N)
    # and its optimal pmf is
    #       pmf(n) = (X_n+1 - X_n)/(X - X_O) \propto (X_n+1 - X_n)
    # Therefore, if we want to mimic this with a simple Geom stop time,
    # we must ensure (X_n+1 - X_n) is decreasing
    # note: let l_n:=log(X_n). Then
    # X_{n+1}-X_n = exp(l_{n+1})-exp(l_n)
    # = exp(l_n)[exp(l_{n+1})exp(-l_n)-1]= exp(l_n)[exp(l_{n+1}-l_n)-1]
    # = exp(l_n)expm1(l_{n+1}-l_n)
    find_offset = function(dta_conv,min_mass){
      # dta_conv=self$dta_adapt[[1]]$dta_trunc
      # find last peak to ensure monotone decreasing
      if(is.null(dta_conv[["cauchy_err"]])){
        cauchy_err=exp(dta_conv$ll_vec[-nrow(dta_conv)])*
          expm1(dta_conv$ll_cauchy_err[-1L])
      }else{cauchy_err=dta_conv$cauchy_err[-1L]}
      cauchy_err[is.nan(cauchy_err) | is.infinite(cauchy_err)]=0
      floor_cerr=pmin(1E-5*max(cauchy_err),.Machine$double.eps)
      sm_ce=pmax(floor_cerr,cauchy_err) # smooth small and <0 blips
      drev=diff(rev(sm_ce))
      ind_origin=length(sm_ce)-which(drev<0)[1L]+1L
      if(is.na(ind_origin)) ind_origin=1L # already monotone decreasing

      # find offset via pmf0 criterion
      # need to be robust to numerical errors that corrupt monotonicity
      pmf_0=exp(dta_conv$ll_vec-dta_conv$ll_vec[nrow(dta_conv)])
      ind_origin=max(ind_origin,which(pmf_0>min_mass)[1L])
      return(ind_origin) # note: returns index not an actual N
    },

    study_joint_convergence = function(
      theta         = self$theta_0,
      min_mass      = 0.8,
      min_eps_speed = if(self$correct_unif) 20.0 else 0.1,
      tol           = sqrt(.Machine$double.eps)
      ){

      # find preliminary offsets and eps_speed, and init storage
      dta_adapt     = self$dta_adapt[[1L]]
      O_trunc       = with(dta_adapt,dta_trunc$N[self$find_offset(dta_trunc,min_mass)])
      O_eps         = with(dta_adapt,dta_eps$N[self$find_offset(dta_eps,min_mass)])
      N_hi_trunc    = max(dta_adapt$dta_trunc$N)
      N_hi_eps      = max(dta_adapt$dta_eps$N)
      N_trunc_vec   = seq.int(O_trunc,N_hi_trunc)
      N_eps_vec     = numeric(length(N_trunc_vec))
      N_eps_vec[1L] = O_eps
      eps_speed     = max(min_eps_speed, (N_hi_eps-O_eps)/length(N_trunc_vec)) # init at default
      ll_vec        = numeric(length(N_trunc_vec))
      ll_vec[1L]    = self$loglik_biased(theta, N_trunc = O_trunc, N_eps = O_eps)
      
      # explore joint convergence, check monotonicity holds (it must except for unif)
      i             = 1L
      N_trunc       = N_trunc_vec[i]
      N_eps         = O_eps
      ind_last_doub = 1L
      ll_true       = max(self$dta_adapt[[1L]]$dta_trunc$ll_vec)
      ll            = ll_vec[1L]
      while(N_trunc<N_hi_trunc && ll+tol<ll_true){
        i       = i + 1L
        N_trunc = N_trunc_vec[i]
        N_eps   = N_eps_vec[i-1L] + eps_speed
        ll      = self$loglik_biased(theta, N_trunc = N_trunc, N_eps = N_eps)
        
        # this loop corrects for the potential non-double-monotonicity of unif
        # by doubling eps_speed until the joint sequence is increasing
        while(ll<ll_vec[i-1L] && ll+tol<ll_true){
          if(self$debug) cat("eps_speed doubled.\n")
          eps_speed     = 2*eps_speed
          ind_last_doub = i
          N_eps         = N_eps_vec[i-1L] + eps_speed
          ll            = self$loglik_biased(theta, N_trunc = N_trunc, N_eps = N_eps)
        }
        
        N_eps_vec[i] = N_eps
        ll_vec[i]    = ll
      }
      
      # build the joint convergence dataframe
      # if we doubled after the first loop step, we need to reconstruct ll_vec
      N_trunc_vec = N_trunc_vec[1L:i]
      if(ind_last_doub<=2L){
        N_eps_vec = N_eps_vec[1L:i]
        ll_vec    = ll_vec[1L:i]
      }else{
        N_eps_vec      = seq(O_eps,by=eps_speed,length.out = i)
        ll_vec_new     = numeric(i)
        ll_vec_new[1L] = ll_vec[1L] # at least the origin is correct
        for(j in seq.int(2L,i)){
          ll_vec_new[j] = self$loglik_biased(
            theta, N_trunc = N_trunc_vec[j], N_eps = N_eps_vec[j]
          )
        }
        ll_vec=ll_vec_new
      }
      dta_joint = data.frame(
        N             = N_trunc_vec,
        N_trunc       = N_trunc_vec,
        N_eps         = N_eps_vec,
        ll_vec        = ll_vec,
        ll_cauchy_err = c(NA_real_,diff(ll_vec))
      )
      
      # compute {X_{n+1}-X_n} (wrongly called "Cauchy" error by me...) seq from {l_n} seq.
      # subtract min(ll_vec) to normalize and avoid underflow in exp()
      dta_joint$cauchy_err=c(
        NA_real_,
        exp(ll_vec[-nrow(dta_joint)]-min(ll_vec))*expm1(dta_joint$ll_cauchy_err[-1L])
      )
      # smooth small and <0 blips
      cauchy_err=dta_joint$cauchy_err[-1L]
      cauchy_err[is.nan(cauchy_err) | is.infinite(cauchy_err)]=0
      floor_cerr=pmin(1E-5*max(cauchy_err),.Machine$double.eps)
      dta_joint$cauchy_err=c(NA_real_,pmax(floor_cerr,cauchy_err))
      
      # store and return eps_speed
      self$dta_adapt[[1L]]$dta_joint = dta_joint
      return(eps_speed)
    },
    
    set_AST = function(
      min_mass      = 0.8,
      slope_alpha   = 0.99,
      min_p_geom    = 0.1,
      max_p_geom    = 0.9
    ){
      # note: joint convergence data is built using a given min_mass so it cannot
      # be reused. it must be recalculated every time we change min_mass.
      eps_speed = self$study_joint_convergence(min_mass = min_mass)
      dta_joint = self$dta_adapt[[1L]]$dta_joint
      
      # find joint offsets
      # NOTE: no need to recompute eps_speed because we set offsets via the
      # the the same N_vec so speed is maintained
      ind_origin=self$find_offset(dta_joint,min_mass)
      O_trunc=dta_joint$N_trunc[ind_origin];O_eps=dta_joint$N_eps[ind_origin]

      # check if we have enough points for regression (>=5)
      if(ind_origin>=nrow(dta_joint)-3L){
        p_geom=max_p_geom # we're already at the very right-end if ind_reg is short
      }else{
        # find decay of geometric stop time by mimicking optimal stopping
        # probability pmf(n) \propto X_n+1 - X_n
        stop_prob=dta_joint$cauchy_err[-1L]
        stop_prob=stop_prob/sum(stop_prob)
        dta_joint$stop_prob=c(stop_prob,0)
        # with(dta_joint[-nrow(dta_joint),],plot(N,stop_prob,log="y",type = "l"))
        # abline(v=O_trunc,lty="dotted")
        # abline(h=stop_prob[ind_origin],lty="dotted")
        # define subset of dta_joint to use for regression
        # then center variables
        ind_reg=seq.int(ind_origin,length(stop_prob)) # from offset to last point
        ind_reg=ind_reg[which(dta_joint$stop_prob[ind_reg]>0)] # need to remove potential 0 introduced by shifting
        dta_lm=dta_joint[ind_reg,];dta_lm$lsp=log(dta_lm$stop_prob)
        dta_lm$lsp=dta_lm$lsp-dta_lm$lsp[1L]; dta_lm$N=dta_lm$N-dta_lm$N[1L]
        fit_lm=lm(lsp~0+N,data=dta_lm) # force 0 intercept to only capture slope
        if(length(ind_reg)<=20L){
          min_slope=-fit_lm$coefficients["N"]
        }else{
          min_slope=-confint(fit_lm,level = slope_alpha)["N",2L]
        }
        # curve({stop_prob[ind_origin]*exp(-min_slope*(x-O_trunc))},
        #       from=O_trunc,to=max(dta_joint$N_trunc)-1L,add=TRUE,lty="dashed")
        p_geom=-unname(expm1(-min_slope))#=1-exp(-min_slope)
        p_geom=max(min_p_geom,min(max_p_geom,p_geom)) # enforce safety limit
      }

      # store results and stop_time
      stop_time = GeometricST$new(
        O_trunc   = O_trunc,
        O_eps     = O_eps,
        prob      = p_geom,
        eps_speed = eps_speed
      )
      self$st_list[[1L]] = stop_time
    }
  )
)
