###############################################################################
# Tune samplers
###############################################################################

#' @export
tune_sampler = function(
  sampler,                                             # an object of class MHSamplerReactNet
  outdir          = getwd(),                           # where should we store the tuned parameters
  n_cores         = parallel::detectCores(),           # number of cores to use
  n_samples       = 10000L,                            # 10000 recommended to detect chains that get stuck due to poor approx.
  varmat_mult_vec = seq(0.5, 3, length.out = 8L),      # seq(0.25,4,length.out = 2L*n_cores)
  sd_thresh_vec   = c(0.1, 0.5*(1:4)),                 # thresholds for std. dev. of log(unbiased likelihood). see Schmon et al. (2020) # c(10^((-5):0),1.5,2)
  min_mass_vec    = c(0.0, 0.01, 0.1, 0.2*(1:4), 0.9), # choose offsets so that X_omega / X >= min_mass # c(0,0.01,0.1,0.2,0.4,0.6,0.8,0.9,0.99,0.999)
  def_min_mass    = max(min_mass_vec)                  # must be set high so that initial run proceeds without issues, even if starting from theta true!
  ){
  
  stopifnot(is(sampler, "MHSamplerReactNet"))
  sampler$debug = TRUE # show more info
  
  #######################################
  # tune stopping times
  #######################################
  
  cat("\nAnalyzing sequence asymptotics and tuning stopping times using a safe strategy\n")
  tic=microbenchmark::get_nanotime()
  sampler$study_convergence()
  gc(verbose=FALSE)
  sampler$set_AST(min_mass = def_min_mass)
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  
  cat("\nInfo about tuned stop-times:\n")
  cat("O_trunc:\n");print(sapply(sampler$st_list,function(l){l$O_trunc}))
  cat("O_eps:\n");print(sapply(sampler$st_list,function(l){l$O_eps}))
  cat("prob:\n");print(sapply(sampler$st_list,function(l){l$prob}))
  cat("eps_speed:\n");print(sapply(sampler$st_list,function(l){l$eps_speed}))
  
  # grow state spaces before continuing to get better time measurements
  cat("\nGrowing spaces first for smoother sampling and more faithful timings\n")
  tic=microbenchmark::get_nanotime()
  for(ss in seq_along(sampler$state_spaces)){
    S_max=sampler$st_list[[ss]]$O_trunc+qgeom(0.9999,sampler$st_list[[ss]]$prob)
    n_spaces = length(sampler$state_spaces[[ss]]) # number of spaces available
    if((missing_spaces <- S_max-(n_spaces-1L))>0L)
      for(m in seq_len(missing_spaces)) sampler$grow_state_spaces(ss)
  }
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  gc(verbose=FALSE)
  
  #######################################
  # preliminary run to improve theta_0 and get variance under posterior
  # use fewer sampler since we are not estimating ESSs here
  #######################################
  
  n_samples_med = as.integer(n_samples/4)
  cat(sprintf("\nShort run of %d samples per chain to get better estimates of mode and variance under posterior\n",n_samples_med))
  tic=microbenchmark::get_nanotime()
  res_run=sampler$run_chains_parallel(n_samples_med, n_cores)
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  res_set=sampler$set_theta_0_from_par_run(res_run)
  cat(sprintf("Estimate of mode with target=%.6f found at:\n",res_set$target))
  print(res_set$theta)
  chainvar=sampler$set_varmat_from_par_run(res_run=res_run) # implicitly mult=1, so this is just variance of the chain
  
  cat("\nAnalyzing truncation and approx. solver asymptotics at new theta_0\n")
  tic=microbenchmark::get_nanotime()
  sampler$study_convergence();gc(verbose=FALSE)
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  
  #######################################
  # estimate sd(loglik) for each configuration
  #######################################
  
  # define list of parameter values to search
  list_pars = list(
    sd_thresh_vec   = sd_thresh_vec,
    varmat_mult_vec = varmat_mult_vec,
    min_mass_vec    = min_mass_vec
  )
  cat("\nWe'll be searching these parameters\n");str(list_pars)
  
  # estimate sd(loglik) and time to compute using iid samples
  cat("\nEstimating sd(loglik) for every configuration.\n")
  tic=microbenchmark::get_nanotime()
  sd_stats = sampler$est_sd_loglik(min_mass_vec=min_mass_vec,n_chains=n_cores)
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  cat("\nEstimated sd(loglik) for each config:\n");print(sd_stats)
  best_pars=sampler$match_sd_thresh_config(sd_stats, sd_thresh_vec)
  cat("\nBest configuration for every threshold on sd(loglik)\n");print(best_pars)
  ind_dup = duplicated(best_pars[,-1]) # check duplicates
  if(any(ind_dup)){
    cat("After removing thresholds with duplicated configurations, we get\n")
    best_pars = best_pars[!ind_dup,]
    print(best_pars)
  }
  
  #######################################
  # select best stoptime configuration with fixed proposal variance = 1*(posterior variance)
  #######################################
  
  tuning_results = varmat_list = vector("list",nrow(best_pars)) # define storage
  for(i in seq_along(tuning_results)){
    tuning_results[[i]]=cbind(best_pars[i,],varmat_mult=1)
    
    # tune stop times to the best config found for this threshold
    cat(sprintf("\nTuning stop times to the best config found for sd(loglik)<=%.1e.\n",
                best_pars$sd_thresh[i]))
    sampler$set_AST(min_mass=best_pars$min_mass[i])
    cat("\nInfo about tuned stop-times:\n")
    cat("O_trunc:\n");print(sapply(sampler$st_list,function(l){l$O_trunc}))
    cat("O_eps:\n");print(sapply(sampler$st_list,function(l){l$O_eps}))
    cat("prob:\n");print(sapply(sampler$st_list,function(l){l$prob}))
    cat("eps_speed:\n");print(sapply(sampler$st_list,function(l){l$eps_speed}))
    
    cat("\nRunning all chains using this configuration to estimate ess/min.\n")
    tic=microbenchmark::get_nanotime()
    res_run=sampler$run_chains_parallel(S=n_samples,n_cores)
    total_sampling_time = ((microbenchmark::get_nanotime() - tic)*1E-9)/60 # minutes
    cat(sprintf("Total sampling time: %.1f minutes\n",total_sampling_time))
    
    # get ess for each chain and parameter and take the mean across all
    ess_mean = mean(sapply(res_run, function(l){
      pmin(coda::effectiveSize(coda::mcmc(l$theta)),mcmcse::ess(l$theta))
    }))
    ess_per_min = ess_mean/total_sampling_time
    cat(sprintf("\nAverage ESS/min: %.2f\n", ess_per_min))
    tuning_results[[i]]$total_sampling_time=total_sampling_time
    tuning_results[[i]]$ess_mean=ess_mean
    tuning_results[[i]]$ess_per_min=ess_per_min
  }
  tuning_results=as.data.frame(do.call("rbind",lapply(tuning_results, unlist)))
  cat("\nFinished! Summary of ESS/min:\n");print(tuning_results)
  
  cat("\nSetting sampler to the best configuration found\n")
  ind_best=which.max(tuning_results$ess_per_min)
  best_tuning_results = tuning_results[ind_best,]
  print(best_tuning_results)
  sampler$set_AST(min_mass=tuning_results$min_mass[ind_best])
  
  #######################################
  # run in parallel to improve estimate of posterior variance
  # this is necessary because the previous run was made using an un-tuned diagonal variance
  #######################################
  
  cat(sprintf("\nShort run of %d samples per chain to get better estimates of mode and variance under posterior\n",n_samples_med))
  tic=microbenchmark::get_nanotime()
  res_run=sampler$run_chains_parallel(n_samples_med, n_cores)
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  res_set=sampler$set_theta_0_from_par_run(res_run)
  cat(sprintf("Estimate of mode with target=%.6f found at:\n",res_set$target))
  print(res_set$theta)
  chainvar=sampler$set_varmat_from_par_run(res_run=res_run) # implicitly mult=1, so this is just variance of the chain
  
  #######################################
  # tune proposal variance for the best stoptime config found
  #######################################
  
  tic=microbenchmark::get_nanotime()
  eff_dta=tryCatch({
    sampler$tune_varmat_mult(
      S               = n_samples,
      varmat_mult_vec = varmat_mult_vec,
      chainvar        = chainvar,
      n_chains        = n_cores
    )
  }, error=function(e){print(e);NULL} # print error info but don't stop execution
  ) 
  cat(sprintf("(Spent %.1f minutes on this task)\n",
              1E-9*(microbenchmark::get_nanotime() - tic)/60))
  if(is.null(eff_dta)){
    cat("\n\nTuning variance multiplier failed. Inspect output above.\n")
  }else{
    best_tuning_results = cbind(best_tuning_results[,1:4], eff_dta)
    cat("\n\nTuning variance multiplier finished. Summary:\n")
    print(best_tuning_results)
  }
  
  #######################################
  # store results
  #######################################
  
  cat("\nStoring results.\n")
  
  # template filename
  res_path = file.path(outdir,paste0(
    sprintf(
      "%s_%s%s",
      sampler$exp_name,
      ifelse(is(sampler, "MHSamplerReactNetRTS"),"reg_ts_",""),
      sampler$gtp_solver
    )
  ))
  if(!dir.exists(res_path)) dir.create(res_path)
  
  # store variance matrix
  write.table(x = sampler$varmat,
              file = file.path(res_path,"varmat.txt"),
              row.names = FALSE,col.names = FALSE)
  
  # store theta_0
  write.table(x = sampler$theta_0,
              file = file.path(res_path,"theta_0.txt"),
              row.names = FALSE,col.names = FALSE)
  
  # store best tuning parameters
  write.table(x = best_tuning_results,
              file = file.path(res_path,"tuning_results.txt"),
              row.names = FALSE,col.names = FALSE)
  
  # store stopping time parameters
  st_dta = do.call(
    rbind,
    lapply(sampler$st_list, function(S){
      data.frame(O_trunc=S$O_trunc, O_eps=S$O_eps, eps_speed=S$eps_speed, prob=S$prob)
    })
  )
  write.table(
    x         = st_dta,
    file      = file.path(res_path, "st_dta.txt"),
    row.names = FALSE
  )
  cat("\nFinished!\n")
}


