##############################################################################
# Metropolis-Hastings sampler for Reaction Network models
# Base class that implements methods dealing with
# - proposal
# - managing collections of statespaces
# - building rate matrices
# - tuning samplers
##############################################################################

MHSamplerReactNet = R6::R6Class(
  classname = "MHSamplerReactNet",
  inherit = MHSampler,
  public = list(
    exp_name          = NULL, # optional: name of the experiment
    dta               = NULL, # list with components s_pre (jump from), s (jump to), dt (time step)
    CTMC              = NULL, # CTMC object representing the reaction network
    n_obs             = NULL, # number of observations
    n_spec            = NULL, # number of species
    n_react           = NULL, # number of reactions
    state_spaces      = NULL, # a list with increasing collections of state spaces. A collection is again a list containing a matrix of states (a state space)
    obs_indices       = NULL, # contains the position of every observation in each statespace
    update_indices    = NULL, # similar to state_spaces, for each state space we have which reaction (if it exists) gets us from each state to another
    up_csv_hash_table = NULL, # "csv" version of the updates stored as hashed environment (fast lookup)
    varmat            = NULL, # variance of proposal
    chmat             = NULL, # Cholesky decomposition of variance matrix of the proposal
    ss_lbound         = NULL, # absolute lower bound for state spaces
    ss_ubound         = NULL, # absolute upper bound for state spaces
    gtp_solver        = NULL, # solver used by get_trans_prob
    correct_unif      = NULL, # logical: should we correct for non-double-monotonicity in sequential uniformization?
    st_list           = NULL, # list of adaptive stop-times, one per observation
    dta_adapt         = NULL, # info about convergence of approx, one per observation
    theta_to_rrs      = NULL, # function: transform theta to reaction rates (defaults to identity)
    verbose           = NULL, # logical, print stoptime? passed to get_trans_prob
    debug             = NULL, # debug mode?
    initialize = function(exp_name = NULL, dta, CTMC, varmat, gtp_solver = "unif",
                          theta_to_rrs = function(x){x},
                          ss_lbound, ss_ubound, verbose = FALSE, debug = FALSE, ...){
      super$initialize(...)
      self$exp_name     = exp_name
      self$dta          = dta
      self$CTMC         = CTMC
      self$n_obs        = nrow(self$dta$t) # remember that dta is a list of matrices!
      self$n_spec       = ncol(self$CTMC$updates)
      self$n_react      = nrow(self$CTMC$updates)
      # if(missing(varmat)) 
      #   varmat          = 0.1*diag(self$theta_0^2)
      self$set_varmat(varmat)
      if(missing(ss_lbound))
        ss_lbound       = rep.int(0L,self$n_spec)
      self$ss_lbound    = ss_lbound
      if(missing(ss_ubound))
        ss_ubound       = rep(Inf,self$n_spec)
      self$ss_ubound    = ss_ubound
      self$theta_to_rrs = theta_to_rrs
      self$verbose      = verbose
      self$debug        = debug
      
      # select matrix exponentiation algorithm
      self$correct_unif   = TRUE
      if(gtp_solver == "wrong_unif"){
        self$correct_unif = FALSE
        gtp_solver        = "unif"
      }
      self$gtp_solver     = gtp_solver
      
      # initialize hash table of updates
      # key: fingerprint of the reaction (csv'd update vector) / value: index of the reaction
      up_csv = apply(self$CTMC$updates, 1L, function(u){ # "csv" the updates
        paste0(u, collapse = ",")})
      up_csv_hash_table = new.env(hash = TRUE, size = self$n_react)
      for(r in seq_len(self$n_react)) up_csv_hash_table[[up_csv[r]]] = r
      self$up_csv_hash_table = up_csv_hash_table
    },
    
    ######################################################################
    # methods for the proposal distribution
    # currently using random walk multivariate normal
    ######################################################################
    
    # define RNG for proposal
    rprop = function(cur){ # multivariate normal
      drop(cur + crossprod(rnorm(self$dim),self$chmat))
    },
    
    # define logdensity of proposal
    ldprop = function(cur, prop){ 0 }, # irrelevant because prop is symmetric
    
    # set varmat and its cholesky decomposition from given matrix
    set_varmat = function(varmat){
      self$varmat = varmat
      # add tolerance to numeric errors by truncating eigenvalues
      eig_dec = eigen(x=varmat,symmetric = TRUE,only.values = FALSE)
      trunc_varmat = tcrossprod(
        eig_dec$vectors * rep(pmax(eig_dec$values,.Machine$double.eps),
                                   each=self$dim),
        eig_dec$vectors)
      self$chmat = chol(trunc_varmat)
    },
    
    # set variance from short parallel run
    set_varmat_from_par_run = function(res_run,n_chains=4L,n_samples=4000L,
                                       varmat_mult=1){
      if(missing(res_run))
        res_run=self$run_chains_parallel(n_samples,n_chains)
      varmat=varmat_mult*var(do.call("rbind",
                                     lapply(res_run,function(l){l$theta})))
      self$set_varmat(varmat)
      invisible(varmat)
    },
    
    ######################################################################
    # methods for computing likelihoods
    ######################################################################
    
    # translate between state of the sampler and the underlying CTMC object
    set_params = function(theta){self$CTMC$react_rates = self$theta_to_rrs(theta)},
    ldprior    = function(...){stop("not implemented")}, # log-density of prior of theta
    loglik     = function(...){stop("not implemented")}, # log likelihood
    ldtarget   = function(theta){                        # get posterior logdensity
      ldp = unname(self$ldprior(theta))     # prior contribution
      if(is.infinite(ldp) && ldp < 0){     
          llik = -Inf                       # avoid computing likelihood when prior=0
      } else {
          llik = unname(self$loglik(theta)) # likelihood contribution
      }
      return(list(target = ldp+llik, loglik = llik, prior = ldp))
    },
    
    ######################################################################
    # methods for building state spaces
    ######################################################################
    
    # for every observed jump, get shortest trajectory between endpoints
    # this initializes the field state_spaces with the element corresponding to
    # K=0 -- i.e., base approximation
    # IDEA: if (s_pre, s) is a valid observation, there must exist a 
    # sequence of updates that joins s_pre and s. Thus, ds:=s-s_pre
    # can be written as a conic combination of the update vectors, 
    # i.e., U^T * x = ds, with x>=0. We get the shortest path by minimizing
    # the L1 norm of x subject to these constraints. This is a linear program!
    # TODO: no support for absolute upper bound
    set_base_state_spaces_shortest_path = function(){
      # solve for the minimal combination of updates that yields given endpoints
      dta=self$dta; n_obs=self$n_obs; n_spec=self$n_spec; n_react=self$n_react
      updates = self$CTMC$updates
      lbound = self$ss_lbound
      ds = matrix(NA_integer_,nrow = n_obs, ncol=n_spec)
      update_comb = matrix(NA_integer_,nrow = n_obs, ncol=n_react)
      trajectories = lapply(seq_len(n_obs), function(i) dta$s_pre[i,,drop=FALSE]) # initialize trajectory with starting point
      t_up = t(self$CTMC$updates)
      for(i in seq_len(n_obs)){
        ds[i,] = dta$s[i,] - dta$s_pre[i,]
        sol_lp = limSolve::linp(E=t_up, F=ds[i,], Cost = rep(1L,n_react))
        update_comb[i,]=sol_lp$X
        stopifnot(!sol_lp$IsError && all.equal( # check it makes sense
          as.vector(t_up %*% update_comb[i,]), ds[i,]))
        
        # complete the trajectory between endpoints with the updates, while
        # ensuring that we do not go outside boundaries
        n_leaps = rep.int(0L, n_react) # steps taken using each update
        cur_st = trajectories[[i]]
        while(!all(n_leaps == update_comb[i,])){
          r=1L # attempt to move with first reaction
          while (n_leaps[r] >= update_comb[i,r] || any(cur_st + updates[r,] < lbound) ) {
            if(r==n_react)
              stop(sprintf("set_base_state_spaces: impossible to not move outside lbound in obs=%d",i))
            else
              r=r+1L
          }
          cur_st = cur_st + updates[r,]
          n_leaps[r] = n_leaps[r] + 1L
          trajectories[[i]] = rbind(trajectories[[i]],cur_st)
        }
        
        # check that last state is in fact the endpoint
        stopifnot(
          all.equal(trajectories[[i]][nrow(trajectories[[i]]),],dta$s[i,]))
      } # end loop observations
      
      # set result as first element of state_spaces (overwrites!)
      self$state_spaces=lapply(trajectories,function(M){list(M)})
    },
    
    # utility: return rows in matrix S that match s
    find_state_in_space = function(s,S){
      which(rowSums(abs(S-rep(s, each=nrow(S)))) < 100*.Machine$double.eps)
    },
    
    # grow a given statespace collection 
    # Adapted from "makeStatespace" function in
    # https://github.com/ageorgou/roulette/
    grow_one_state_space_col = function(ss_col,all_ways=TRUE){
      n_spec=self$n_spec; n_react=self$n_react
      lbound = self$ss_lbound; ubound = self$ss_ubound
      if(all_ways){ # move in all directions
        I_spec = matrix(as.integer(diag(self$n_spec)),self$n_spec) # force integerness
        updates = rbind(I_spec,-1L*I_spec)
      }else{ # move using the reactions of the CTMC
        updates = self$CTMC$updates
      }
      n_spaces=length(self$state_spaces[[ss_col]])
      current_space = self$state_spaces[[ss_col]][[n_spaces]] # extract most recent
      new_states = do.call("rbind",
        lapply(seq_len(nrow(updates)),function(r){ # loop updates
          # apply the r-th update to all states in current_space
          current_space + rep(updates[r,],each=nrow(current_space))
        })) # rbind all the states we just created
      # check which states don't exceed limits
      ind_keep = apply(new_states,1L,function(s){all(s>=lbound & s<=ubound)})
      new_space = unique( # filter, append, remove duplicates
        rbind(current_space, new_states[ind_keep,,drop=FALSE]))
      self$state_spaces[[ss_col]][[n_spaces+1L]] = new_space # add to collection
    },
    
    # this fn creates the update index for the most recent state space in a
    # given state-space collection
    # for a given state space, the "update index" is a triplet list
    # (s0,s1,r) that tells us reaction r takes us *directly* from s0 to s1,
    # if there is such a reaction
    populate_update_index = function(ss_col){
      up_csv=self$up_csv; up_csv_hash_table=self$up_csv_hash_table
      # extract most recent state space in the collection
      n_spaces=length(self$state_spaces[[ss_col]])
      st_space = self$state_spaces[[ss_col]][[n_spaces]]
      
      # need to compute the pairwise jump matrix and then "csv" it
      # start by getting all jump changes for all species
      # the "-" effectively "transposes" the result
      l_mat = lapply(seq_len(self$n_spec), function(i){
        -outer(st_space[,i],st_space[,i],"-")
      })
      # intersperse l_mat with commas in order to csv them together coordinate-wise
      comma_list = vector("list", 2L*self$n_spec-1L)
      for (i in seq_along(comma_list)) {
        if(i%%2L==0L)
          comma_list[[i]]=","
        else
          comma_list[[i]]=l_mat[[i%/%2L+1L]]
      }
      # csv!
      ds_csv = matrix(do.call("paste0", comma_list),nrow=nrow(st_space))
      # TODO: this is still pretty slow
      # note: the approach above is MUCH faster than this misleadingly short code
      # ds_csv =apply(st_space, 1L, function(m){ # "csv" the jumps
      #       apply(st_space, 1L, function(n){paste0(m-n, collapse = ",")})
      #     })
      # find the csv'd jumps in the hash table to check if there's such a reaction
      # hashing is >14 times faster than using match!
      match_vec = unlist(mget(x=ds_csv,envir = up_csv_hash_table, ifnotfound = 0L),
                         use.names = FALSE)
      ind_nz = which(match_vec>0L)-1L # 0-based
      ind_nz_i = ind_nz %% nrow(st_space); ind_nz_j = ind_nz %/% nrow(st_space)
      # build the triplet list
      # note: dims is necessary otherwise constructor might chop trailing
      # rows or cols that only have zeroes
      if(packageVersion("Matrix") >= "1.3.0")
        new_update_index = Matrix::sparseMatrix(
          i = ind_nz_i,j = ind_nz_j, dims=rep.int(nrow(st_space),2L),
          index1 = FALSE,x = match_vec[ind_nz+1L], repr = "T") # keep triplet format
      else
        new_update_index = Matrix::sparseMatrix(
          i = ind_nz_i,j = ind_nz_j, dims=rep.int(nrow(st_space),2L),
          index1 = FALSE,x = match_vec[ind_nz+1L], giveCsparse=FALSE) # keep triplet format
      # append update index to the ss_col-th collection
      self$update_indices[[ss_col]][[n_spaces]] = new_update_index
    },
    
    ######################################################################
    # methods for building rate matrix
    ######################################################################

    # build rate matrix: computed once per likelihood evaluation :(
    # note: the Q's produced are not conservative because we don't add a
    # "phantom" state to collect the missing flux generated by truncating
    # This does not affect the transition probabilities estimated!
    Q_fun = function(ss_col,K){
      # ss_col: index of the collection of state spaces
      # K: approximation order (0 based)
      up_ind_mat = self$update_indices[[ss_col]][[K+1L]] # 1-based R indexing
      st_space =  self$state_spaces[[ss_col]][[K+1L]] # 1-based R indexing
      
      # find nonzero transitions and compute corresponding rates
      nz_vals = mapply(function(x,i){
        self$CTMC$get_propensity(x,st_space[i,])
      }, up_ind_mat@x,up_ind_mat@i+1L) # Matrix uses 0-based indexing
      
      # get diagonal = negative sum of the rates of all possible reactions
      dQ = -colSums(apply(st_space, 1L, self$CTMC$get_all_propensities))
      
      # build sparse matrix from triplets
      # need as.numeric for when nz_vals is empty list, since
      # unlist(list())==NULL (bad) but as.numeric(list())==numeric(0) (good)
      n = nrow(st_space); seq_n=seq.int(0L,n-1L)
      Q = Matrix::sparseMatrix(
        i=c(up_ind_mat@i,seq_n), j=c(up_ind_mat@j,seq_n), dims = rep.int(n,2L),
        x = c(as.numeric(nz_vals),dQ), index1 = FALSE) # use 0-based index
      
      if(any(Matrix::rowSums(Q)>1E-7)){ # rowSums<=0 except for rounding error
        warning(
          sprintf("Q matrix inconsistent for ss_col=%d and K=%d, since max(rowSums(Q)) = %g.\n",
                  ss_col,K,max(Matrix::rowSums(Q))))
      }
      return(Q)
    },
    
    # build rate matrix including a phantom state so that resulting Q is 
    # conservative.
    Q_fun_complete = function(ss_col,K){
      # ss_col: index of the collection of state spaces
      # K: approximation order (0 based)
      up_ind_mat = self$update_indices[[ss_col]][[K+1L]] # 1-based R indexing
      st_space =  self$state_spaces[[ss_col]][[K+1L]] # 1-based R indexing
      
      # find nonzero transitions and compute corresponding rates
      nz_vals = mapply(function(x,i){
        self$CTMC$get_propensity(x,st_space[i,])
      }, up_ind_mat@x,up_ind_mat@i+1L) # Matrix uses 0-based indexing
      
      # get diagonal = negative sum of the rates of all possible reactions
      dQ = -colSums(apply(st_space, 1L, self$CTMC$get_all_propensities))
      
      # build sparse matrix from triplets
      # need as.numeric for when nz_vals is empty list, since
      # unlist(list())==NULL (bad) but as.numeric(list())==numeric(0) (good)
      n = nrow(st_space); seq_n=seq.int(0L,n-1L)
      Q = Matrix::sparseMatrix(
        i=c(up_ind_mat@i,seq_n), j=c(up_ind_mat@j,seq_n), dims = rep.int(n+1L,2L),
        x = c(as.numeric(nz_vals),dQ), index1 = FALSE) # use 0-based index
      Q[,n+1L] = -Matrix::rowSums(Q) # add flux to absorbing state
      if(any(Matrix::rowSums(Q)>100*.Machine$double.eps)){ # rowSums<=0 except for rounding error
        warning(
          sprintf("Q matrix inconsistent for ss_col=%d and K=%d, since max(rowSums(Q)) = %g.\n",
                  ss_col,K,max(Matrix::rowSums(Q))))
      }
      return(Q)
    },
    
    ######################################################################
    # methods for tuning parameters
    ######################################################################
    
    # tune m in var(proposal) = m*var(chain)
    tune_varmat_mult = function(varmat_mult_vec,n_chains,chainvar,
                                S=2000L, min_acc=0.1, thrsh_ess_ratio=0.15){
      res = parallel::mclapply(X=varmat_mult_vec,FUN = function(mult){
        sink(file = self$log_file, append = TRUE)
        self$set_varmat(chainvar*mult)
        tic=microbenchmark::get_nanotime()
        res=self$run_chain(S)
        toc=1E-9*(microbenchmark::get_nanotime()-tic)/60 # minutes
        acc_ratio = 1-res$stats$n_reject/S
        ind=which.max(res$ld_stats[,"target"])
        # we take min of coda and mcmcse because the latter is overly optimistic
        # when the chain is sticky (common when sampler has not been tuned properly yet)
        ess = pmin(
          coda::effectiveSize(coda::mcmc(res$theta)),
          mcmcse::ess(res$theta))
        mean_ess = mean(ess)
        sink()
        list(varmat_mult=mult,acc_ratio=acc_ratio,theta=res$theta[ind,],
             target=unname(res$ld_stats[ind,"target"]),mean_ess=mean_ess,
             ess_per_min=mean_ess/toc)
      }, mc.cores = n_chains)
      ind_err=which(sapply(res, is, class2 = "try-error")) # any errors?
      if(length(ind_err)>0L) stop(res[[ind_err[1L]]]) # break on first error
      
      # find variance maximizing ess/min, as long as it yields acc_ratio>=min_acc,
      # and it yields at most thrsh_ess_ratio
      # too high ess suggests the estimate is bad and that the pseudomarginal method failed
      eff_dta=as.data.frame(
        do.call("rbind",lapply(res, function(l){l$theta=NULL;unlist(l)})))
      ind_cand = ((eff_dta$acc_ratio>=min_acc) & 
                    (eff_dta$mean_ess/S < thrsh_ess_ratio))
      if(!any(ind_cand)){
        warn="No varmat_mult achieves minimum quality thresholds. ESS estimates are not credible so setting to 0.\n"
        cat(warn);warning(warn)
        best_config=eff_dta[which.min(min_acc-eff_dta$acc_ratio),]
        best_config$ess_per_min = 0
      }else{
        candidates=eff_dta[ind_cand,]
        best_config=candidates[which.max(candidates$ess_per_min),]
      }
      self$set_varmat(chainvar*best_config$varmat_mult)
      cat("Selected the following configuration:\n");print(best_config)
      invisible(best_config)
    },
    
    # estimate std dev of the residual log likelihood and time to compute for
    # a given min_mass_vec.
    # Note: no need to call study_convergence() since that will not change
    # with min_mass (but dta_joint will!)
    est_sd_loglik = function(min_mass_vec, R=200L, n_chains, alpha=0.99){
      res_list=parallel::mclapply(
        X=seq_along(min_mass_vec),FUN=function(r){
          self$set_AST(min_mass=min_mass_vec[r])
          tic = microbenchmark::get_nanotime()
          ll_vec = replicate(R,self$loglik(self$theta_0),simplify = FALSE)
          time_spent = (microbenchmark::get_nanotime()-tic)*1E-9
          ll_vec = unlist(ll_vec)
          if(any(is.infinite(ll_vec))){
            sd_loglik=Inf
          }else{
            # get traditional sd and one implied by quantiles
            sd_loglik = sd(ll_vec);mean_ll=mean(ll_vec)
            quants=quantile(ll_vec,c((1-alpha)/2,(1+alpha)/2))
            norm_q=qnorm((1+alpha)/2)
            sd_loglik_quants=max(quants[2L] - mean_ll, mean_ll-quants[1L])/norm_q
            sd_loglik=max(sd_loglik,sd_loglik_quants)
          }
          c("min_mass"=min_mass_vec[r],"sd_loglik"=sd_loglik, 
            "time_per_loglik"=time_spent/R)
        }, mc.cores = n_chains, mc.preschedule = FALSE)
      ind_err=which(sapply(res_list, is, class2 = "try-error")) # any errors?
      if(length(ind_err)>0L){ # print errors, break on first
        print(res_list[ind_err]);stop(res_list[[ind_err[1L]]])
      }
      sd_stats=as.data.frame(do.call("rbind",res_list))
      sd_stats$sd_loglik[sd_stats$sd_loglik < 100*.Machine$double.eps] = NA_real_ # filter 0-ish variances
      return(sd_stats)
    },
    # find the best config for each threshold for sd(loglik)
    match_sd_thresh_config = function(par_grid, sd_thresh_vec){
      if(all(is.na(par_grid$sd_loglik)))
        stop("All configurations yield invalid estimates of sd(loglik).\n")
      par_grid = par_grid[!is.na(par_grid$sd_loglik),]
      best_pars = data.frame()
      for(i in seq_along(sd_thresh_vec)){ # loop thresholds
        sd_thresh = sd_thresh_vec[i]
        candidates = par_grid[par_grid$sd_loglik<=sd_thresh,] # find candidate configurations
        if(nrow(candidates) == 0L){
          warn=sprintf(
            "No configuration achieves sd(loglik)<=%.1e. Using the closest match.\n",
            sd_thresh)
          cat(warn);warning(warn)
          best_pars = rbind(best_pars,par_grid[which.min(par_grid$sd_loglik-sd_thresh),])
        } else if(nrow(candidates) == 1L){
          best_pars = rbind(best_pars,candidates)
        } else{ # select the fastest
          best_pars = rbind(best_pars,candidates[which.min(candidates$time_per_loglik),])
        }
      }
      best_pars=cbind(sd_thresh=sd_thresh_vec,best_pars)
      return(best_pars)
    }
  )
)
