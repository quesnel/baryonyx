library(rbaryonyx)

f <- function(x) {
  r <- rbaryonyx::optimize_01lp_problem(file_path="bibd1n.lp",
           limit = 8000,
           theta = 0.5,
           delta = 0.01,
           constraint_order = 1,
           kappa_min = x,
           kappa_step = 1e-10,
           kappa_max = 2.0,
           alpha = 2.0,
           w = 60.0,
           time_limit = 10.0,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

result <- optimize(f, lower=1, upper=1, maximum=FALSE, tol=0.00001)

result <- optim(par=0.3, f, lower=0, upper=1, method="Brent")

###
###

# optimizer returns: 3910.000000 (kappa: 0.116151 0.000000 5.000000
#                                 delta: 0.000705 theta: 0.784961)
# optimizer returns: 3860.000000 (kappa: 0.118276 0.000001 5.000000
#                                 delta: 0.000615 theta: 0.784961)

# theta 7.849608e-01
# delta 6.145989e-04
# kappa_min 1.182756e-01
# kappa_max 8.983414e-07

library(rgenoud)
library(rilp)
library(parallel)

optim_gen_lp <- function(x) {
  r <- rilp::optimize_01lp_problem(
           file_path = "/home/gquesnel/devel/lp/verger_10_10-opt.lp",
           limit = 4000,
           theta = x[1],
           delta = x[2],
           constraint_order = 1,
           kappa_min = x[3],
           kappa_step = x[4],
           kappa_max = 5.0,
           alpha = 1.0,
           w = 60,
           time_limit = 30,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

d = matrix(c(0.5, 0.0005,  0.1, 1e-7,
             0.9, 0.001,  0.15, 9e-7),
             nrow=4, ncol=2)

s = c(0.8, 0.003226, 0.118947, 2e-7)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile="debug.txt")

claw1 <- genoud(optim_gen_lp, nvars=4,
                Domains=d,
                starting.values=s,
                cluster=cl,
                max=FALSE, pop.size=10)
