library(rbaryonyx)
library(sensitivity)

factors = c(
  "theta",
  "delta",
  "constraint_order",
  "kappa_min",
  "kappa_step",
  "alpha",
  "w",
  "norm",
  "float_type"
  )

bounds = data.frame(
  min=c(
    0,     # theta
    0,     # delta
    0,     # constraint_order
    0,     # kappa_min
    1e-16, # kappa_step
    0.0,   # alpha
    50,    # w
    0,     # norm
    0.0,   # init_policy_random
    0.0,   # init_random
    0
    ),    # float_type
max=c(
    1,     # theta
    0,     # delta
    4,     # constraint_order
    0.1,   # kappa_min
    1e-1,  # kappa_step
    2.0,   # alpha
    500,   # w
    4,     # norm
    1.0,   # init_policy_random
    1.0,   # init_random
    2))    # float_type

rownames(bounds) <- factors

morrisDesign <- morris(model = NULL,
                factors = factors,
                r = 10,
                design=list(type="oat", levels=10, grid.jump=5),
                binf = bounds$min,
                bsup = bounds$max,
                scale=TRUE)

# morrisDesign$X
#
# 20 * (8 + 1)
#

solve_lp <- function(x, file_path, limit=10000, time_limit=10, seed=123456789, thread=1) {
  r <- rbaryonyx::solve_01lp_problem(file_path = file_path,
                   limit = limit,
                   theta = x["theta"],
                   delta = x["delta"],
                   constraint_order = x["constraint_order"],
                   kappa_min = x["kappa_min"],
                   kappa_step = x["kappa_step"],
                   alpha = x["alpha"],
                   w = x["w"],
                   time_limit = time_limit,
                   seed = seed,
                   thread = thread,
                   norm = x["norm"],
                   init_policy_random = x["init_policy"],
                   init_random = x["init_random"],
                   float_type = x["float_type"])

  return ifelse(r$solution_found, r$value, 10000000)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="verger_5_5.lp", thread=1, limit=10000, time_limit=10, seed=123456789)

morrisDesign$Y <- r
mu <- apply(morrisDesign$X,2,mean)
mu.star <- apply(morrisDesign$X, 2, function(x) mean(abs(x)))
sigma <- apply(morrisDesign$ee, 2, sd)

apply(morrisDesign$X, 2, function(v) plot(factor(v), r))
