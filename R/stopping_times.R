#######################################
# Stopping time
#######################################

StoppingTime = R6::R6Class(
  classname = "StoppingTime",
  public = list(
    set_params = function(...){
      alist=list(...)
      for(p in names(alist)) self[[p]] = alist[[p]]
    },
    rtime = function(...){stop("not implemented")}, # RNG
    pmf = function(...){stop("not implemented")},
    lpmf = function(x){log(self$pmf(x))}, # default, usually specialized methods do better
    cdf = function(...){stop("not implemented")},
    surv = function(...){stop("not implemented")} # 1-cdf
  )
)


#######################################
# Geometric
#######################################

GeometricST = R6::R6Class(
  classname = "GeometricST",
  inherit = StoppingTime,
  public = list(
    prob=NULL,
    O_trunc=NULL,
    O_eps=NULL,
    eps_speed=NULL,
    initialize = function(prob,O_trunc,O_eps,eps_speed){
      self$prob=prob
      self$O_trunc=O_trunc
      self$O_eps=O_eps
      self$eps_speed=eps_speed
    },
    rtime = function(n=1){
      N=rgeom(n = n,prob = self$prob)
      return(list(N=N,N_trunc=self$O_trunc+N,N_eps=self$O_eps+N*self$eps_speed,
                  pmf=self$pmf(N), lpmf=self$lpmf(N)))
    },
    pmf = function(x){
      dgeom(x=x,prob = self$prob)
    },
    lpmf = function(x){
      dgeom(x=x,prob = self$prob,log = TRUE)
    },
    cdf = function(q){
      pgeom(q=q,prob=self$prob)
    },
    surv = function(q){
      pgeom(q=q,prob=self$prob,lower.tail = FALSE)
    },
    mean = function(){
      p = self$prob
      return((1-p)/p)
    }
  )
)

# # test
# geom_st = GeometricST$new(0.5)
# geom_st$rtime()
# geom_st$surv(0:5)
