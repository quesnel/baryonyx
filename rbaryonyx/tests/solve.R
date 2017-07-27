library(rbaryonyx)
library(sensitivity)

factors=c("theta", "delta", "constraint_order", "kappa_min",
                   "kappa_step", "kappa_max", "alpha", "w")
bounds = data.frame(min=c(0,   0, 0, 0.0,    0,  1.0, 1.0, 50),
                    max=c(1, 0.1, 4, 0.25, 0.1,  1.0, 1.0, 50))
rownames(bounds) <- factors

morrisDesign <- morris(model = NULL,
                factors = factors,
                r = 100,
                design=list(type="oat", levels=20, grid.jump=5),
                binf = bounds$min,
                bsup = bounds$max,
                scale=TRUE)

# morrisDesign$X
#
# 20 * (8 + 1)
#

solve_lp <- function(x, file_path, thread=1, limit=10000, time_limit=10) {
  theta <- x["theta"]
  delta <- x["delta"]
  constraint_order <- x["constraint_order"]
  kappa_min <- x["kappa_min"]
  kappa_step <- x["kappa_step"] 
  kappa_max <- x["kappa_max"] 
  alpha <- x["alpha"]
  w <- x["w"]

  r <- rbaryonyx::optimize_01lp_problem(file_path = file_path,
           limit = limit,
           theta = theta,
           delta = delta,
           constraint_order = constraint_order,
           kappa_min = kappa_min,
           kappa_step = kappa_step,
           kappa_max = kappa_max,
           alpha = alpha,
           w = w,
           time_limit = time_limit,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="/home/gquesnel/devel/lp/verger_10_10-opt.lp")
mr = matrix(unlist(r), ncol=2, byrow = TRUE)
# mr1 <- mr[,-3]

tell(morrisDesign, mr1)
plot(morrisDesign)


##
## Try to search optimal parameteres for the verger_10_10-opt.lp 
## 
##

library(rbaryonyx)
library(sensitivity)

factors=c("theta", "delta", "kappa_min",
                   "kappa_step")

bounds = data.frame(min=c(0.3, 0.0001, 0.10, 1e-10),
                    max=c(  1, 0.1,    0.25, 1e-3))

rownames(bounds) <- factors

morrisDesign <- morris(model = NULL,
                factors = factors,
                r = 100,
                design=list(type="oat", levels=20, grid.jump=5),
                binf = bounds$min,
                bsup = bounds$max,
                scale=TRUE)

solve_lp <- function(x, file_path, thread=1, limit=10000, time_limit=10) {
  theta <- x["theta"]
  delta <- x["delta"]
  kappa_min <- x["kappa_min"]
  kappa_step <- x["kappa_step"] 

  r <- rbaryonyx::optimize_01lp_problem(file_path = file_path,
           limit = limit,
           theta = theta,
           delta = delta,
           constraint_order = 1,
           kappa_min = kappa_min,
           kappa_step = kappa_step,
           kappa_max = 1.0,
           alpha = 1.0,
           w = 60,
           time_limit = time_limit,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="/home/gquesnel/devel/lp/verger_10_10-opt.lp")
mr = matrix(unlist(r), ncol=2, byrow = TRUE)
# mr1 <- mr[,-3]

tell(morrisDesign, mr1)
plot(morrisDesign)
