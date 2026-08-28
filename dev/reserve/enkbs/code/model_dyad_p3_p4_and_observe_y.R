## acir reserve file
## Origin: aci/R/benchmark_models.R:88-176
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Whole aci model_dyad body: the p3 preset is jiang2026enkbs; the p4 preset and observe='y' branch are cross-filed (see the POINTER files in paper-extremes/ and extensions/).
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The whole aci 0.0.30 model_dyad body, retained so the p3 (jiang2026enkbs) and p4 (moser2026extremes) presets and the observe=='y' non-CGNS branch (package extension) can be re-applied verbatim. The acir mainline keeps only the p1 / observe=='x' arms.

model_dyad <- function(variant = c("p1", "p3", "p4"),
                       observe = c("x", "y"), params = NULL) {
  variant <- match.arg(variant); observe <- match.arg(observe)
  defaults <- switch(variant,
    p1 = list(d_x = 0.5, gamma = 2,   f_x = 0.5, s_x = 0.5, d_y = 0.5, f_y = 1,   s_y = 1),
    p3 = list(d_x = 0.5, gamma = 2,   f_x = 1,   s_x = 0.5, d_y = 0.5, f_y = 0.8, s_y = 1),
    p4 = list(d_x = 0.8, gamma = 1.2, f_x = 1,   s_x = 0.5, d_y = 0.8, f_y = 1,   s_y = 2))
  p <- .complete_scalar_params(params, defaults, "dyad",
                               c("d_x", "gamma", "s_x", "d_y", "s_y"))
  if (observe == "x") {
    # The observed-x dyad has a package batch realiser.  Its coefficient
    # closures share a locked environment so the attached realiser descriptor
    # cannot silently outlive mutation of captured constructor parameters.
    coefficient_env <- list2env(list(p = p), parent = baseenv())
    lockEnvironment(coefficient_env, bindings = TRUE)
    coefficient_functions <- list(
      Lx = function(t, x) matrix(p$gamma * x, 1, 1),
      fx = function(t, x) -p$d_x * x + p$f_x,
      Ly = function(t, x) matrix(-p$d_y, 1, 1),
      fy = function(t, x) -p$gamma * x^2 + p$f_y,
      Sx1 = function(t, x) matrix(p$s_x, 1, 1),
      Sx2 = function(t, x) matrix(0, 1, 1),
      Sy1 = function(t, x) matrix(0, 1, 1),
      Sy2 = function(t, x) matrix(p$s_y, 1, 1)
    )
    coefficient_functions <- lapply(coefficient_functions, function(fun) {
      environment(fun) <- coefficient_env
      fun
    })
    m <- cgns_model(
      Lx = coefficient_functions$Lx,
      fx = coefficient_functions$fx,
      Ly = coefficient_functions$Ly,
      fy = coefficient_functions$fy,
      Sx1 = coefficient_functions$Sx1,
      Sx2 = coefficient_functions$Sx2,
      Sy1 = coefficient_functions$Sy1,
      Sy2 = coefficient_functions$Sy2,
      k = 1, l = 1, name = sprintf("dyad[%s] y->x", variant))
  } else {
    m <- stochastic_model(
      f = function(t, x, y) -p$d_y * x - p$gamma * y^2 + p$f_y,
      g = function(t, x, y) (-p$d_x + p$gamma * x) * y + p$f_x,
      Sx = function(t, x) matrix(p$s_y, 1, 1),
      Sy = function(t, x, y) matrix(p$s_x, 1, 1),
      k = 1, l = 1, vectorized_members = TRUE,
      name = sprintf("dyad[%s] x->y (non-injective u^2)", variant))
  }
  m$meta$energy_conserving <- TRUE
  m$meta$params <- p
  m$meta$vars <- if (observe == "x")
    list(observed = "x", hidden = "y") else
    list(observed = "y", hidden = "x")
  m$meta$provenance <- switch(variant,
    p1 = paste("andreou2026aci Sections 3.1 and SI.4.1;",
               "ACI_code-main/dyad_interaction_model.m"),
    p3 = "jiang2026enkbs Section 3.2 bidirectional nonlinear-dyad EnKBS causal-inference benchmark",
    p4 = "moser2026extremes equations (4.1)-(4.2)")
  m$meta$source_status <- if (observe == "y" && variant == "p3") {
    "paper + MATLAB checked (published EnKBS dyad experiment)"
  } else if (observe == "y") {
    paste("package extension (source equations; reverse x -> y partition",
          "is not a supplied paper/MATLAB inference benchmark)")
  } else if (variant == "p1") {
    "paper + MATLAB checked"
  } else if (variant == "p3") {
    paste("paper checked; the published EnKBS MATLAB studies this direction",
          "with the ensemble engine")
  } else {
    "paper checked; no corresponding MATLAB supplied"
  }
  m$meta$partition_caveat <- if (observe == "y")
    paste("The reverse x -> y partition is not CGNS because the hidden x",
          "enters through x^2; this constructor is for ensemble methods.",
          if (variant == "p3")
            "jiang2026enkbs Section 3.2 studies this direction with EnKBS." else "") else NULL
  m$meta$anti_damping_threshold <- p$d_x / p$gamma
  # The MATLAB reference starts each component at its uncoupled forced
  # equilibrium (F_x / d_x, F_y / d_y).  In particular, the default
  # andreou2026aci configuration starts at (1, 2), not (1, 0).
  m$meta$ic_default <- if (observe == "x")
    list(x0 = p$f_x / p$d_x, y0 = p$f_y / p$d_y) else
    list(x0 = p$f_y / p$d_y, y0 = p$f_x / p$d_x)
  if (observe == "x")
    m <- .attach_cgns_realizer(
      m, "dyad_observed_x_v1", list(params = p)
    )
  m
}
