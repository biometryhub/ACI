# The grading matrix -----------------------------------------------------------
#
# This file exists because of a defect this package found in itself, and it is
# the mechanism that makes that class of defect impossible rather than merely
# remembered.
#
# The defect: a term can be publicly reachable, execute on every run, be
# reported as covered by line coverage AND be graded by an oracle, and still
# be checked against nothing, because the fixture that grades it happens to
# multiply it by zero. That is what happened to the scalar noise
# cross-covariance, and it is the failure mode of every oracle whose scope is
# narrower than the surface it appears to grade.
#
# Remembering to check is not a control. So this file MEASURES: for every
# component of the schema and every consumer of it, it builds the components
# each grading scenario actually uses, records whether the term is non-zero
# there, and fails if a reachable cell is neither graded nor declared as a
# known gap with a reason.
#
# Adding a term, a consumer, or a fixture that silently annihilates something
# therefore breaks the build. That is the point: the register cannot be
# satisfied by intending to check.

# ---- The surface being tracked ----------------------------------------------

.aci_gm_terms <- c(
  "L_x", "f_x", "L_y", "f_y",
  "S_xoS_x", "S_yoS_y", "S_yoS_x"
)

# Consumers, and whether each exists for a given schema variant. A consumer
# that does not exist cannot have an ungraded term; a consumer that exists and
# is reachable must have every term it touches accounted for.
.aci_gm_consumers <- list(
  scalar = c("filter", "smoother", "metric", "online_smoother", "cir"),
  vector = c("filter", "smoother", "metric", "online_smoother", "cir")
)

# ---- Known gaps, each with a reason ------------------------------------------
#
# A cell may appear here only with a stated reason. This list records what the
# package does NOT grade, written down where a reader will find it rather than
# discovered later.
# Empty, and that is a measured statement rather than an aspiration. This list
# held two entries when the register was first built. The online smoother and
# the causal influence range were graded only by a fixture whose noise
# cross-covariance was zero. The register is what found them; the
# `cir_cross` fixture is what closed them.
.aci_gm_known_gaps <- list()

# ---- Grading scenarios -------------------------------------------------------
#
# Each entry is a components list that an oracle test ACTUALLY grades, paired
# with the consumers that oracle compares against a reference. These mirror the
# oracle test files; if they ever drift apart the union below changes and the
# assertions notice.

.aci_gm_fixture <- function(name) {
  system.file("extdata", name, package = "aciR")
}

.aci_gm_scenarios <- function() {
  dyad_x <- read.csv(.aci_gm_fixture("dyad_signal_x.csv"), header = FALSE)$V2
  p <- list(d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1,
            sigma_x = 0.5, sigma_y = 1)

  cross_x <- read.csv(.aci_gm_fixture("cross_signal_x.csv"), header = FALSE)$V2
  cross <- list(
    L_x = 2 * cross_x, f_x = 0.5 - 0.5 * cross_x, L_y = -0.5,
    f_y = 1 - 2 * cross_x^2,
    S_xoS_x = 0.6^2 + 0.3^2, S_yoS_y = 0.5^2 + 0.8^2,
    S_yoS_x = 0.5 * 0.6 + 0.8 * 0.3, S_xoS_y = 0.5 * 0.6 + 0.8 * 0.3
  )

  pp <- read.csv(.aci_gm_fixture("predprey_signal.csv"))
  pp_params <- list(alpha = 0.4, beta = 0.1, gamma = 1.1, delta = 0.4,
                    sigma_x = 0.3, sigma_y = 0.3)

  mv_sig <- read.csv(.aci_gm_fixture("mv_signal.csv"))
  s <- matrix(
    c(0.60, 0.10, 0.25, 0.05, 0.20, 0.50, 0.10, 0.30,
      0.15, 0.05, 0.70, 0.10, 0.05, 0.10, 0.15, 0.55), 4L, 4L
  )
  s_x <- s[1:2, ]
  s_y <- s[3:4, ]
  n_mv <- nrow(mv_sig)
  mv <- list(
    L_x = array(1, c(2L, 2L, n_mv)), f_x = matrix(1, 2L, n_mv),
    L_y = array(1, c(2L, 2L, n_mv)), f_y = matrix(1, 2L, n_mv),
    S_xoS_x = s_x %*% t(s_x), S_yoS_y = s_y %*% t(s_y),
    S_yoS_x = s_y %*% t(s_x)
  )

  list(
    list(id = "dyad", variant = "scalar",
         comp = aci_dyad_components(dyad_x, p),
         consumers = c("filter", "smoother", "metric")),
    list(id = "cross", variant = "scalar", comp = cross,
         consumers = c("filter", "smoother", "metric")),
    list(id = "predprey_a", variant = "scalar",
         comp = aci_predprey_components(pp$prey, pp_params, "predator_to_prey"),
         consumers = c("filter", "smoother", "metric")),
    list(id = "predprey_b", variant = "scalar",
         comp = aci_predprey_components(pp$predator, pp_params,
                                        "prey_to_predator"),
         consumers = c("filter", "smoother", "metric")),
    list(id = "cir", variant = "scalar",
         comp = aci_dyad_components(dyad_x[seq_len(2001L)], p),
         consumers = c("online_smoother", "cir")),
    list(id = "cir_cross", variant = "scalar",
         comp = local({
           cc <- cross
           cc$L_x <- cross$L_x[seq_len(2001L)]
           cc$f_x <- cross$f_x[seq_len(2001L)]
           cc$f_y <- cross$f_y[seq_len(2001L)]
           cc
         }),
         consumers = c("online_smoother", "cir")),
    list(id = "mv_online", variant = "vector", comp = mv,
         consumers = c("online_smoother", "cir")),
    list(id = "mv", variant = "vector", comp = mv,
         consumers = c("filter", "smoother", "metric")),
    list(id = "mv_collapse", variant = "vector",
         comp = local({
           dx <- dyad_x[seq_len(1500L)]
           nn <- length(dx)
           list(L_x = array(2 * dx, c(1L, 1L, nn)),
                f_x = matrix(0.5 - 0.5 * dx, 1L, nn),
                L_y = matrix(-0.5, 1L, 1L),
                f_y = matrix(1 - 2 * dx^2, 1L, nn),
                S_xoS_x = matrix(0.25, 1L, 1L),
                S_yoS_y = matrix(1, 1L, 1L),
                S_yoS_x = matrix(0, 1L, 1L))
         }),
         consumers = c("online_smoother", "cir")),
    list(id = "enso", variant = "vector",
         comp = local({
           e <- read.csv(.aci_gm_fixture("enso_signal.csv"))
           aci_enso_components(e$T_C, e$T_E, e$I, time = e$t)
         }),
         consumers = c("filter", "smoother", "metric"))
  )
}

# A term is exercised only if it is present AND not identically zero. "Present
# but zero" is precisely the state that reads as covered and grades nothing.
.aci_gm_exercised <- function(comp, term) {
  value <- comp[[term]]
  !is.null(value) && any(abs(value) > 0)
}

.aci_gm_build <- function() {
  scenarios <- .aci_gm_scenarios()
  rows <- list()
  for (variant in names(.aci_gm_consumers)) {
    for (consumer in .aci_gm_consumers[[variant]]) {
      for (term in .aci_gm_terms) {
        graded <- any(vapply(scenarios, function(s) {
          identical(s$variant, variant) && consumer %in% s$consumers &&
            .aci_gm_exercised(s$comp, term)
        }, logical(1)))
        rows[[length(rows) + 1L]] <- data.frame(
          variant = variant, consumer = consumer, term = term,
          graded = graded, stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

.aci_gm_is_known_gap <- function(variant, consumer, term) {
  any(vapply(.aci_gm_known_gaps, function(g) {
    identical(g$variant, variant) && identical(g$consumer, consumer) &&
      identical(g$term, term)
  }, logical(1)))
}

# ---- The assertions ----------------------------------------------------------

test_that("every reachable term is graded, or declared with a reason", {
  matrix_df <- .aci_gm_build()

  ungraded <- matrix_df[!matrix_df$graded, ]
  undeclared <- ungraded[!mapply(
    .aci_gm_is_known_gap, ungraded$variant, ungraded$consumer, ungraded$term
  ), ]

  if (nrow(undeclared) > 0L) {
    message(
      "Ungraded and undeclared cells:\n",
      paste(
        sprintf(
          "  %s / %s / %s", undeclared$variant, undeclared$consumer,
          undeclared$term
        ),
        collapse = "\n"
      )
    )
  }
  expect_identical(nrow(undeclared), 0L)
})

test_that("every declared gap is really a gap", {
  # A stale exemption is worse than none: it licenses a term that has since
  # been graded, and would go on licensing it if the grade were later lost.
  matrix_df <- .aci_gm_build()

  # Asserted before the loop, so this test still says something when the
  # exemption list is empty. A test whose only assertions live inside a loop
  # over an empty list is a test that skips, and a gate that skips is not a
  # gate, which is the same standard the oracle tests are held to.
  expect_identical(sum(!matrix_df$graded), length(.aci_gm_known_gaps))

  for (gap in .aci_gm_known_gaps) {
    row <- matrix_df[
      matrix_df$variant == gap$variant &
        matrix_df$consumer == gap$consumer &
        matrix_df$term == gap$term,
    ]
    expect_identical(nrow(row), 1L)
    expect_false(
      row$graded,
      info = sprintf(
        "%s/%s/%s is declared a gap but is now graded; remove the exemption.",
        gap$variant, gap$consumer, gap$term
      )
    )
    expect_gt(nchar(gap$reason), 40L)
  }
})

test_that("the register covers the schema it claims to", {
  # If a component is added to the schema without being added here, the
  # register would silently stop describing the surface. Hold it to the
  # documented schema instead of to itself.
  documented <- c(
    "L_x", "f_x", "L_y", "f_y",
    "S_xoS_x", "S_yoS_y", "S_yoS_x", "S_xoS_y"
  )
  # S_xoS_y is the transpose of S_yoS_x and is graded with it, so it is
  # tracked through that entry rather than separately.
  expect_setequal(c(.aci_gm_terms, "S_xoS_y"), documented)

  matrix_df <- .aci_gm_build()
  expect_identical(
    nrow(matrix_df),
    length(.aci_gm_terms) *
      sum(vapply(.aci_gm_consumers, length, integer(1)))
  )
  # And the register must not be trivially satisfiable: at least one cell has
  # to be genuinely graded, or every assertion above passes vacuously.
  expect_gt(sum(matrix_df$graded), 20L)
})

# The sibling class, declared time-varying and read as constant ----------------
#
# The register above catches a term that is present but ANNIHILATED. It cannot
# catch a term that is present, non-zero, and silently read only at its first
# step, which produces a plausible result computed from an input the caller
# did not supply.
#
# That defect was real. The vector validator inverted only the first slice of a
# time-varying observation-noise covariance, so a covariance ramping from 0.5
# to 0.9 gave a result identical, to the last bit, to a constant 0.5. No oracle
# caught it: the fixtures all hold that covariance constant. No contract caught
# it: the array was accepted. Line coverage was total throughout.
#
# The general property is stronger than non-degeneracy and subsumes it:
#
#   a coefficient is genuinely consumed only if PERTURBING IT CHANGES THE
#   ANSWER.
#
# So each term that may vary in time is perturbed at a LATE step (late enough
# that a reader of the first step alone would not notice), and the output must
# move. A term whose perturbation changes nothing is either ignored or dead.

.aci_gm_perturb_at <- function(n) as.integer(round(0.75 * n))

.aci_gm_scalar_system <- function(n = 300L) {
  dt <- 1e-3
  x <- 1 + 0.3 * sin(5 * seq(0, by = dt, length.out = n))
  comp <- list(
    L_x = 0.8 + 0.2 * x, f_x = 0.4 - 0.3 * x, L_y = rep(-0.6, n),
    f_y = 0.5 - 0.2 * x^2,
    S_xoS_x = 0.5, S_yoS_y = 0.9, S_yoS_x = 0.2, S_xoS_y = 0.2
  )
  list(x = x, comp = comp, dt = dt, n = n)
}

.aci_gm_vector_system <- function(n = 300L) {
  dt <- 1e-3
  tt <- seq(0, by = dt, length.out = n)
  x <- rbind(1 + 0.3 * sin(5 * tt), 0.5 + 0.2 * cos(7 * tt))
  arr <- function(m) array(rep(m, n), c(2L, 2L, n))
  comp <- list(
    L_x = arr(matrix(c(0.8, 0.1, 0.2, 0.7), 2L, 2L)),
    f_x = matrix(0.2, 2L, n),
    L_y = arr(matrix(c(-1.1, 0.2, 0.3, -0.9), 2L, 2L)),
    f_y = matrix(0.3, 2L, n),
    S_xoS_x = arr(matrix(c(0.7, 0.1, 0.1, 0.8), 2L, 2L)),
    S_yoS_y = arr(matrix(c(0.9, 0.2, 0.2, 1.0), 2L, 2L)),
    S_yoS_x = arr(matrix(c(0.2, 0.05, 0.05, 0.15), 2L, 2L))
  )
  list(x = x, comp = comp, dt = dt, n = n)
}

.aci_gm_filter_out <- function(sys) {
  filt <- aci_filter(
    sys$x, sys$comp, sys$dt,
    mu0 = if (is.matrix(sys$x)) c(0.4, 0.1) else 0.4,
    R0 = if (is.matrix(sys$x)) diag(c(0.3, 0.3)) else 0.3
  )
  c(as.numeric(filt$mean), as.numeric(filt$cov))
}

test_that("perturbing a late step of any time-varying term moves the result", {
  for (sys in list(.aci_gm_scalar_system(), .aci_gm_vector_system())) {
    base <- .aci_gm_filter_out(sys)
    j <- .aci_gm_perturb_at(sys$n)

    for (term in .aci_gm_terms) {
      value <- sys$comp[[term]]
      d <- dim(value)
      varying <- (is.null(d) && length(value) == sys$n) ||
        (length(d) == 2L && ncol(value) == sys$n) ||
        (length(d) == 3L && d[3L] == sys$n)
      if (!varying) {
        next
      }

      bumped <- sys
      if (length(d) == 3L) {
        # Scale the whole slice, preserving symmetry so a covariance stays a
        # covariance and the perturbation tests consumption rather than
        # admissibility.
        bumped$comp[[term]][, , j] <- value[, , j] * 1.5
      } else if (length(d) == 2L) {
        bumped$comp[[term]][, j] <- value[, j] + 0.5
      } else {
        bumped$comp[[term]][j] <- value[j] + 0.5
      }

      moved <- max(abs(.aci_gm_filter_out(bumped) - base))
      expect_gt(moved, 0)
    }
  }
})

test_that("a time-varying observation covariance is not read as constant", {
  # The specific defect, pinned so it cannot return. Both systems agree at the
  # FIRST step and differ later; a filter that reads only the first step gives
  # bit-identical answers, which is what this once did.
  n <- 300L
  dt <- 1e-3
  tt <- seq(0, by = dt, length.out = n)
  x <- rbind(1 + 0.3 * sin(5 * tt), 0.5 + 0.2 * cos(7 * tt))
  build <- function(varying) {
    s_xx <- array(0, c(2L, 2L, n))
    for (j in seq_len(n)) {
      first <- if (varying) 0.5 + 0.4 * (j - 1L) / (n - 1L) else 0.5
      s_xx[, , j] <- diag(c(first, 0.8))
    }
    list(
      L_x = diag(2L), f_x = c(0, 0), L_y = -diag(2L), f_y = c(0, 0),
      S_xoS_x = s_xx, S_yoS_y = diag(2L), S_yoS_x = matrix(0, 2L, 2L)
    )
  }
  expect_identical(build(TRUE)$S_xoS_x[, , 1L], build(FALSE)$S_xoS_x[, , 1L])

  a <- aci_filter(x, build(TRUE), dt, mu0 = c(0, 0), R0 = diag(2L))
  b <- aci_filter(x, build(FALSE), dt, mu0 = c(0, 0), R0 = diag(2L))
  expect_gt(max(abs(a$cov - b$cov)), 1e-6)
})

test_that("a covariance that loses positivity later is rejected, not ignored", {
  # The same first-slice-only reading would have accepted a covariance that is
  # positive definite at step one and singular afterwards.
  n <- 200L
  x <- matrix(0, 2L, n)
  s_xx <- array(0, c(2L, 2L, n))
  for (j in seq_len(n)) {
    s_xx[, , j] <- diag(c(if (j < n %/% 2L) 0.5 else -0.2, 0.8))
  }
  comp <- list(
    L_x = diag(2L), f_x = c(0, 0), L_y = -diag(2L), f_y = c(0, 0),
    S_xoS_x = s_xx, S_yoS_y = diag(2L), S_yoS_x = matrix(0, 2L, 2L)
  )
  expect_error(
    aci_filter(x, comp, 1e-3, mu0 = c(0, 0), R0 = diag(2L)),
    "positive definite at every"
  )
})
