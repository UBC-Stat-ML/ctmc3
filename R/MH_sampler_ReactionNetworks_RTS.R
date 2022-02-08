##############################################################################
# Metropolis-Hastings sampler for Reaction Network models
# Special case for Regular Time Series (RTS):
#   - state_spaces has a unique increasing collection of state spaces
#   - this collection is used to compute trans-prob for all observations
##############################################################################

MHSamplerReactNetRTS = R6::R6Class(
  classname = "MHSamplerReactNetRTS",
  inherit = MHSamplerReactNet,
  public = list(
    initialize = function(dta=dta,...) {
      vdt = as.numeric(var(dta$dt))
      if(vdt > .Machine$double.eps)
        stop("Not a regularly sampled time series.")
      super$initialize(dta=dta,...)
      if(!self$debug) self$init_state_spaces() # initialize state spaces
    },
    
    ######################################################################
    # methods for computing likelihoods
    ######################################################################
    
    # compute (biased) estimate of all transition probabilites
    get_all_trans_probs_est = function(theta,N_trunc,N_eps){
      # N_trunc/N_eps: truncation/solver order of approximation
      # if(any(theta<0)) return(-Inf)# WRONG: THIS SHOULD BE HANDLED BY THE PRIOR (TODO!)
      self$set_params(theta) # translate theta to model's parameters, and set it
      if(self$verbose) cat(sprintf("loglik-bias: N_trunc=%d, N_eps=%.1f\n",
                                   N_trunc,N_eps))
      
      # check if N_trunc-th state space is available, grow if needed
      n_spaces = length(self$state_spaces[[1L]])
      if((missing_spaces<-N_trunc - (n_spaces-1L))>0L)
        for(m in seq_len(missing_spaces)) self$grow_state_spaces() # grow if needed
      
      Q_N = self$Q_fun(1L, N_trunc) # get Q matrix at N_trunc-th order of approximation
      
      # get positions of the observation in the state space
      n_N = self$obs_indices[[N_trunc+1L]]$s_pre
      m_N = self$obs_indices[[N_trunc+1L]]$s
      
      # compute probabilities using solver
      p_N = get_trans_prob(n = n_N, m=m_N, dt=self$dta$dt[1L],Q = Q_N,ms = 1L,
                           eps = 10^(-N_eps), solver=self$gtp_solver,
                           verbose = self$verbose)
      return(p_N)
    },
    
    # log of (biased) likelihood estimate
    loglik_biased = function(theta,N_trunc,N_eps){
      sum(log(self$get_all_trans_probs_est(theta,N_trunc,N_eps)))
    },
    
    # get log of unbiased estimate of likelihood
    # we de-bias the *product of all terms* without explicitly computing products
    # Want: log(Z) where Z = prod(p_O)+(prod(pN1)-prod(pN))/pmf(N)
    # Factorizing:
    # Z = prod(p_O)[1+(prod(pN1)-prod(pN))/(pmf(N)prod(p_O))]
    # Define S_n:=sum(log(p_n)). Then
    # log(Z) = S_O+log1p[(prod(pN1)-prod(pN))/(pmf(N)prod(p_O))]
    # Also
    # (prod(pN1)-prod(pN))/(pmf(N)prod(p_O))
    # = [exp(S_N1)-exp(S_N)]/[pmf(N)exp(S_O)]
    # = exp(S_N-S_O)[exp(S_N1-S_N)-1]/pmf(N)
    # = exp(S_N-S_O-log(pmf(N)))expm1(S_N1-S_N) =: w
    # Finally, log(Z) = S_O + log1p(w)
    # We need to be able to handle cases where S = -Inf
    # When pN1 =0 => pO=pN=0 (by monotonicity) => ans = -Inf
    # If pN1>0 but pO=pN=0, we have Z = pN1/prob => logZ = SN1-log(pmf(N))
    # If pN1>pN>0 but pO=0, we have Z = (pN1-pN)/prob = (exp(SN1)-exp(SN))/prob
    # = exp(SN)(exp(SN1-SN)-1)/prob = exp(SN-log(prob))expm1(SN1-SN)
    # Hence, logZ = SN - log(prob) + log(expm1(SN1-SN))
    # Note: we use max(0,) to deal with numerical inaccuracies
    loglik = function(theta){
      # if(any(theta<0)) return(-Inf)# WRONG: THIS SHOULD BE HANDLED BY THE PRIOR (TODO!)
      stop_time=self$st_list[[1L]]
      rtime = stop_time$rtime(); N=rtime$N
      N_trunc=rtime$N_trunc; N_eps=rtime$N_eps; lprob = rtime$lpmf
      if(self$verbose) cat(sprintf("loglik: N=%d, N_trunc=%d, N_eps=%.1f\n",
                                   N, N_trunc,N_eps))
      # note: p_X's are vectors
      p_N  = self$get_all_trans_probs_est(theta,N_trunc,N_eps)
      p_N1 = self$get_all_trans_probs_est(theta,N_trunc+1L,N_eps+stop_time$eps_speed)
      S_N  = sum(log(p_N))
      S_N1 = sum(log(p_N1))
      if(is.infinite(S_N1) && S_N1 < 0) return(-Inf)       # means that all 3 are 0 (by domination)
      if(is.infinite(S_N)  && S_N  < 0) return(S_N1-lprob) # Z= 0 + (pN1-0)/prob = pN1/prob
      if(N==0L){
        p_O = p_N
        S_O = S_N
      }else{
        p_O = self$get_all_trans_probs_est(theta,stop_time$O_trunc,stop_time$O_eps)
        S_O = sum(log(p_O))
        if(is.infinite(S_O) && S_O < 0) 
          return(S_N - lprob + log(expm1(max(0,S_N1-S_N))))
      }
      # compute log(Z) when all terms are finite
      w    = exp(max(0, S_N - S_O) - lprob) * expm1(max(0, S_N1 - S_N))
      logZ = S_O + log1p(w)
      return(logZ)
    },
    
    ######################################################################
    # methods for building state spaces
    ######################################################################
    
    # initialize state spaces and associated data structures
    init_state_spaces = function(){
      self$set_base_state_spaces_shortest_path() # shortest path state space
      self$collapse_state_spaces()
      self$obs_indices=list();self$update_indices=vector("list",1L) # reset obs & update indices
      self$build_obs_index()
      self$populate_update_index(1L)
    },
    
    # merges all observations' initial state spaces into an unique state space
    collapse_state_spaces = function(){
      collapsed_ss = unique(
        do.call("rbind",lapply(self$state_spaces, function(l){l[[1L]]})))
      self$state_spaces=list(list(collapsed_ss))
    },
    
    # find the position of every observation in the state spaces
    build_obs_index = function(){
      n_spaces = length(self$state_spaces[[1L]])
      st_space = self$state_spaces[[1L]][[n_spaces]] # extract collapsed state space
      obs_ind_s_pre = apply(self$dta$s_pre,1L,function(s){
        self$find_state_in_space(s,st_space)})
      obs_ind_s = apply(self$dta$s,1L,function(s){
        self$find_state_in_space(s,st_space)})
      self$obs_indices[[n_spaces]] = list(s_pre=obs_ind_s_pre,s=obs_ind_s)
    },
    
    # dispatch a method for growing state spaces
    grow_state_spaces = function(...){
      self$grow_one_state_space_col(1L)
      self$build_obs_index()
      self$populate_update_index(1L)
    },

    ######################################################################
    # methods for stopping times
    ######################################################################

    study_convergence = function(theta=self$theta_0,
                                 max_err_thresh = sqrt(.Machine$double.eps),
                                 min_ll_thresh = -700,
                                 N_lo_eps = if(self$gtp_solver=="skeletoid") -10L else 0L){
      cat("\nStudying convergence of the likelihood\n")

      # set N_eps high and study truncation Cauchy error
      N_hi_eps=ceiling(-log10(.Machine$double.eps))
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
