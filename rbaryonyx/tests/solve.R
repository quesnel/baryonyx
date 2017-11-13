library(rbaryonyx)
library(sensitivity)

factors = c("theta", "delta", "constraint_order", "kappa_min", "kappa_step",
  "kappa_max", "alpha", "w", "norm", "pushing_k_factor",
  "pushing_objective_amplifier", "pushes_limit", "pushing_iteration_limit",
  "float_type")

bounds = data.frame(
  min=c(
    0,     # theta
    0,     # delta
    0,     # constraint_order
    0,     # kappa_min
    1e-16, # kappa_step
    1.0,   # kappa_max
    0.0,   # alpha
    50,    # w
    0,     # norm
    0.1,   # pushing_k_factor
    1.0,   # pushing_objective_amplifier
    10,    # pushes_limit
    20,    # pushing_iteration_limit
    0),    # float_type
max=c(
    1,     # theta
    0,     # delta
    4,     # constraint_order
    0.1,   # kappa_min
    1e-1,  # kappa_step
    1.0,   # kappa_max
    2.0,   # alpha
    500,   # w
    4,     # norm
    1,     # pushing_k_factor
    10.0,  # pushing_objective_amplifier
    100,   # pushes_limit
    200,   # pushing_iteration_limit
    2))    # float_type

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

solve_lp <- function(x, file_path, limit=10000, time_limit=10, seed=123456789, thread=1) {
  r <- rbaryonyx::solve_01lp_problem(file_path = file_path,
                   limit = limit,
                   theta = x["theta"],
                   delta = x["delta"],
                   constraint_order = x["constraint_order"],
                   kappa_min = x["kappa_min"],
                   kappa_step = x["kappa_step"],
                   kappa_max = x["kappa_max"],
                   alpha = x["alpha"],
                   w = x["w"],
                   time_limit = time_limit,
                   seed = seed,
                   thread = thread,
                   norm = x["norm"],
                   pushing_k_factor = x["pushing_k_factor"],
                   pushing_objective_amplifier = x["pushing_objective_amplifier,"],
                   pushes_limit = x["pushes_limit"],
                   pushing_iteration_limit = x["pushing_iteration_limit"],
                   float_type = x["float_type"])

  return(r)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="flat30-7.lp", thread=1, limit=10000, time_limit=10, seed=123456789)

tell(morrisDesign, r)
plot(morrisDesign)
