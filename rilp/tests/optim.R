library(rilp)

f <- function(x) {
  r <- rilp::optimize_01lp_problem(file_path="bibd1n.lp",
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

