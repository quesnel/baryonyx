library(rgenoud)
library(rbaryonyx)
library(parallel)

optim_gen_lp <- function(x) {
  r <- rbaryonyx::optimize_01lp_problem(
           file_path = "verger_10_10-opt.lp",
           limit = 10000,
           theta = x[1],
           delta = x[2],
           constraint_order = 0,
           kappa_min = x[3],
           kappa_step = x[4],
           kappa_max = 1.0,
           alpha = 1.0,
           w = 5000,
           time_limit = 120,
           seed = 123654785,
           thread = 4,
           norm = 4,
           pushing_k_factor = 1,
           pushing_objective_amplifier = 10,
           pushes_limit = 20,
           pushing_iteration_limit = 50,
           init_policy = 0,
           init_policy_random = 0.5,
           init_random = 0.5,
           float_type = 1)
           verbose = FALSE)

  return(r)
}

d = matrix(c(0.25, 0.0005, 0.0, 1e-10,
             0.90, 0.001,  0.2, 1e-6),
             nrow=4, ncol=2)

s = c(0.5, 0.003226, 0.11, 2e-7)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile="debug.txt")

claw1 <- genoud(optim_gen_lp, nvars=4,
                Domains=d,
                starting.values=s,
                cluster=cl,
                max=FALSE, pop.size=10)
