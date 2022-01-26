###############################################################################
# utility to get P_n(X(dt)=m)
###############################################################################

get_trans_prob = function(n,m,dt,Q,ms,eps,verbose=FALSE,solver="skeletoid"){
  # ms: used to convert from state to matrix index
  # solver: "skeletoid", "unif"
  if(verbose) cat(sprintf("\tgtp(%s): n=%d, m=%d, dt=%.2f\n",solver,n,m,dt))
  n_adj = n-ms+1L; m_adj = m-ms+1L # translate state to matrix index
  len_n = length(n); seq_n = seq_along(n)
  v = Matrix::sparseMatrix(i=seq_n,j=n_adj,dims = c(len_n,nrow(Q)),x = rep(1,len_n))
  ind_res = unname(cbind(seq_n,m_adj)) # extracts only the requested jump probabilities

  # pick solver and proceed
  switch (solver,
    skeletoid = {
      pske::skeletoid_vtexpm(Q=Q,t_pow=dt,v=v,eps=eps,verbose=verbose)[ind_res]
    },
    unif = {
      pske::unif_vtexpm(Q=Q,t_pow=dt,v=v,eps=eps,verbose=verbose)[ind_res]
    },
    { stop(sprintf("Unknown solver '%s'.", solver)) } # default case
  )
}
