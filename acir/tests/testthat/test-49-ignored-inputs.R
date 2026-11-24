# Inputs that used to be ignored are refused ---------------------------------

quietly_ii <- function(expr) withCallingHandlers(
  expr, warning = function(w) invokeRestart("muffleWarning"))

ii_closure_model <- function() aci_model(
  Lx = function(t, x) matrix(c(1, 0.5), 1, 2),
  fx = function(t, x) -0.5 * x,
  Ly = function(t, x) diag(-0.5, 2),
  fy = function(t, x) c(0, 0),
  Sx1 = function(t, x) matrix(0.5, 1, 1),
  Sy2 = function(t, x) diag(2), k = 1, l = 2)

ii_setup <- function() {
  dyad <- aci_dyad_model()
  closure <- ii_closure_model()
  list(dyad = dyad,
       dyad_obs = as_obs(simulate(dyad, seed = 2, t_end = 1, dt = 0.01)),
       dyad_init = list(mean = 0, cov = matrix(1)),
       closure = closure,
       closure_obs = as_obs(simulate(closure, seed = 3, t_end = 1, dt = 0.01)),
       closure_init = list(mean = c(0, 0), cov = diag(2)))
}


test_that("aci(table = ) refuses stepper, nsub and regularize", {
  s <- ii_setup()
  for (route in c("dyad", "closure")) {
    m <- s[[route]]; ob <- s[[paste0(route, "_obs")]]
    ini <- s[[paste0(route, "_init")]]
    tb <- lag_table(m, ob, init = ini)
    path <- tempfile(fileext = ".rds"); saveRDS(tb, path)
    tb_rds <- readRDS(path); unlink(path)
    for (table in list(tb, tb_rds)) {
      e <- expect_error(aci(m, ob, table = table, init = ini, nsub = 4),
                        class = "aci_error_dims")
      expect_match(conditionMessage(e), "lag_table()", fixed = TRUE)
      expect_match(conditionMessage(e), "lag table: nsub.", fixed = TRUE)
      expect_error(aci(m, ob, table = table, init = ini, nsub = "x"),
                   class = "aci_error_dims")
      expect_error(aci(m, ob, table = table, init = ini, stepper = "implicit"),
                   class = "aci_error_dims")
      expect_error(aci(m, ob, table = table, init = ini, stepper = "explicit"),
                   class = "aci_error_dims")
      expect_error(aci(m, ob, table = table, init = ini, nsub = 1L),
                   class = "aci_error_dims")
      expect_error(aci(m, ob, table = table, init = ini, regularize = "floor"),
                   class = "aci_error_dims")
      e <- expect_error(aci(m, ob, table = table, init = ini,
                            stepper = "explicit", regularize = "none"),
                        class = "aci_error_dims")
      expect_match(conditionMessage(e), "lag table: stepper, regularize.",
                   fixed = TRUE)
    }
    expect_identical(aci(m, ob, table = tb, init = ini)$aci, lt_diag(tb),
                     info = route)
  }
})


test_that("unknown entries in init are refused on every verb", {
  s <- ii_setup()
  for (route in c("dyad", "closure")) {
    m <- s[[route]]; ob <- s[[paste0(route, "_obs")]]
    ini <- s[[paste0(route, "_init")]]
    bad <- list(m = 5 + ini$mean, cov = ini$cov)
    e <- expect_error(aci_filter(m, ob, init = bad), class = "aci_error_dims")
    expect_match(conditionMessage(e), "unknown entries: m.", fixed = TRUE)
    expect_match(conditionMessage(e), "accepted entries are mean and cov",
                 fixed = TRUE)
    expect_error(aci_filter(m, ob, init = bad, stepper = "implicit"),
                 class = "aci_error_dims")
    expect_error(aci(m, ob, init = bad), class = "aci_error_dims")
    expect_error(aci(m, ob, init = bad, keep = "table"),
                 class = "aci_error_dims")
    expect_error(aci_smoother(m, ob, init = bad), class = "aci_error_dims")
    expect_error(aci_online(m, ob, lag = 3, init = bad),
                 class = "aci_error_dims")
    expect_error(lag_table(m, ob, init = bad), class = "aci_error_dims")

    f <- aci_filter(m, ob, init = ini)
    tb <- lag_table(m, ob, init = ini)
    expect_error(aci_smoother(m, ob, filter = f, init = bad),
                 class = "aci_error_dims")
    expect_error(aci_online(m, ob, lag = 3, filter = f, init = bad),
                 class = "aci_error_dims")
    expect_error(lag_table(m, ob, filter = f, init = bad),
                 class = "aci_error_dims")
    expect_error(aci(m, ob, table = tb, init = bad), class = "aci_error_dims")

    e <- expect_error(aci_filter(m, ob, init = unname(ini)),
                      class = "aci_error_dims")
    expect_match(conditionMessage(e), "(unnamed), (unnamed)", fixed = TRUE)
    expect_error(aci_filter(m, ob, init = c(ini, covariance = 3)),
                 class = "aci_error_dims")
  }
})


test_that("unknown init entries are refused on a library matrix model", {
  me <- aci_enso_model(hidden = c("u", "hW", "tau"))
  oe <- as_obs(simulate(me, seed = 12, t_end = 0.1, dt = 0.005))
  e <- expect_error(
    aci_filter(me, oe, init = list(mu = c(1, 1, 1),
                                   cov = diag(c(0.0032, 0.0008, 0.28)))),
    class = "aci_error_dims")
  expect_match(conditionMessage(e), "unknown entries: mu.", fixed = TRUE)
})


test_that("a valid or stored prior gives the same numbers as before", {
  s <- ii_setup()
  for (route in c("dyad", "closure")) {
    m <- s[[route]]; ob <- s[[paste0(route, "_obs")]]
    ini <- s[[paste0(route, "_init")]]
    f <- aci_filter(m, ob, init = ini)
    stored <- f$meta$init
    expect_identical(names(stored), c("mean", "cov", "source"))
    path <- tempfile(fileext = ".rds"); saveRDS(stored, path)
    stored_rds <- readRDS(path); unlink(path)
    for (prior in list(stored, stored_rds)) {
      expect_identical(aci_filter(m, ob, init = prior)$mean, f$mean)
      expect_identical(aci(m, ob, init = prior)$aci, aci(m, ob, init = ini)$aci)
      expect_identical(lt_diag(lag_table(m, ob, filter = f, init = prior)),
                       lt_diag(lag_table(m, ob, init = ini)))
    }
    a <- aci(m, ob, init = ini, keep = "none")
    expect_no_error(quietly_ii(aci_range(a)))
    expect_identical(
      quietly_ii(aci_filter(m, ob, init = list()))$mean,
      quietly_ii(aci_filter(m, ob))$mean)
  }
})
