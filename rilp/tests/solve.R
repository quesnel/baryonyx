library(rilp)
library(sensitivity)

factors=c("theta", "delta", "constraint_order", "kappa_min",
                   "kappa_step", "kappa_max", "alpha", "w")
bounds = data.frame(min=c(0, 0, 0, 0.0,   0,  1.0, 0.0, 0),
                    max=c(1, 1, 4, 1.0, 0.1, 10.0, 2.0, 1000))
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

f <- function(x, thread=1, limit=1000000, time_limit=1) {
  theta <- x["theta"]
  delta <- x["delta"]
  constraint_order <- x["constraint_order"]
  kappa_min <- x["kappa_min"]
  kappa_step <- x["kappa_step"] 
  kappa_max <- x["kappa_max"] 
  alpha <- x["alpha"]
  w <- x["w"]

  r <- rilp::solve_01lp_problem(file_path="n-queens-problem-20.lp",
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
           verbose = TRUE)

  return(r)
}

r = apply(morrisDesign$X, 1, f)
mr = matrix(unlist(r), ncol=2, byrow = TRUE)
mr1 <- mr[,-3]

tell(morrisDesign, mr1)
plot(morrisDesign)
