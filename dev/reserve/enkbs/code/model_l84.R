## acir reserve file
## Origin: aci/R/benchmark_models.R:321-388
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Default variant 'p3' is jiang2026enkbs; the 'fbcir' arm is cross-filed (see fbcir/code/model_l84_fbcir_variant.POINTER.md). No ACI_code arm exists.
## Verbatim copy from the aci 0.0.30 sources; not modified.

#' Stochastic Lorenz-84; observed (y, z), hidden x -> CGNS
#'
#' @param params Optional complete parameter list overriding the preset.
#' @param variant jiang2026enkbs constant-forcing or andreou2026cir/FBCIR
#'   seasonal preset.
#' @param target In the FBCIR preset, observed target (`"y"` or `"z"`); the
#'   other observed component is recorded as the conditioning channel.
#'
#' @references
#' Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
#' identifying forward and backward causal influence ranges using assimilative
#' causal inference. arXiv:2510.21889v2, 4 August 2026.
#' \doi{10.48550/arXiv.2510.21889}
#'
#' Jiang, Z., Andreou, M., Reich, S. and Chen, N. (2026). A continuous-time
#' ensemble Kalman-Bucy smoother for causal inference and model discovery.
#' arXiv:2604.25157. \doi{10.48550/arXiv.2604.25157}
#' @examples
#' model_l84(variant = "fbcir")
#'
#' @export
model_l84 <- function(params = NULL, variant = c("p3", "fbcir"),
                      target = c("y", "z")) {
  variant <- match.arg(variant)
  target <- match.arg(target)
  p <- params %||% switch(variant,
    p3 = list(a = 0.25, b = 4, F0 = 8, F_amp = 0,
              F_period = 73, g = 1, sig = 0.1),
    fbcir = list(a = 0.25, b = 4, F0 = 8, F_amp = 3,
                 F_period = 73, g = 1, sig = 0.2))
  # Accept the historical custom-parameter spelling Ff while representing the
  # FBCIR seasonal forcing explicitly as F0 + F_amp cos(2 pi t/F_period).
  forcing <- function(t) {
    if (!is.null(p$Ff)) {
      if (is.function(p$Ff)) p$Ff(t) else p$Ff
    } else {
      p$F0 + p$F_amp * cos(2 * pi * t / p$F_period)
    }
  }
  m <- cgns_model(
    Lx = function(t, x) matrix(c(-p$b * x[2] + x[1], p$b * x[1] + x[2]), 2, 1),
    fx = function(t, x) c(p$g - x[1], -x[2]),
    Ly = function(t, x) matrix(-p$a, 1, 1),
    fy = function(t, x) p$a * forcing(t) - x[1]^2 - x[2]^2,
    Sx1 = function(t, x) diag(p$sig, 2),
    Sy2 = function(t, x) matrix(p$sig, 1, 1),
    k = 2, l = 1, name = sprintf("L84[%s] (y,z | x)", variant))
  m$meta$params <- p
  m$meta$variant <- variant
  is_seasonal <- is.null(p$Ff) && !is.null(p$F_amp) && p$F_amp != 0
  m$meta$forcing <- list(type = if (is_seasonal)
    "seasonal_cosine" else "constant_or_custom",
    F0 = p$F0 %||% NA_real_, F_amp = p$F_amp %||% NA_real_,
    F_period = p$F_period %||% NA_real_)
  m$meta$provenance <- if (variant == "fbcir")
    "FBCIR_code-main/lorenz84.m" else "jiang2026enkbs model discovery benchmark"
  m$meta$vars <- list(observed = c("y", "z"), hidden = "x")
  if (variant == "fbcir") {
    m$meta$target_obs_idx <- match(target, c("y", "z"))
    m$meta$conditioning_obs_idx <- setdiff(1:2, m$meta$target_obs_idx)
    m$meta$causal_link <- sprintf("x -> %s | %s", target,
                                  setdiff(c("y", "z"), target))
  }
  # Physical MATLAB ordering is (x, y, z) = (1, 0, 1); this constructor's
  # observed vector is (y, z) and hidden scalar is x.
  m$meta$ic_default <- list(x0 = c(y = 0, z = 1), y0 = c(x = 1))
  m
}
