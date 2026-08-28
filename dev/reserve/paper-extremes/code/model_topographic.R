## acir reserve file
## Origin: aci/R/benchmark_models.R:657-750
## Source package: aci 0.0.30, git tree 97f6b124
## Category: paper-extremes
## Intended release: after the three MATLAB-backed families
## Reason: moser2026extremes eqs. (4.5)-(4.13); spectral barotropic model, not the FBCIR layered one.
## Verbatim copy from the aci 0.0.30 sources; not modified.

#' Spectral barotropic topographic flow (moser2026extremes eqs. 4.5-4.13); V
#' observed.
#' @param noise_convention Complex-noise component convention.
#' @param params Optional complete parameter list.
#'
#' @references
#' Moser, C., Chen, N. and Andreou, M. (2026). Mechanisms and pathways of
#' extreme events in partially-observed stochastic dynamical systems.
#' arXiv:2605.22692. \doi{10.48550/arXiv.2605.22692}
#' @examples
#' model_topographic()
#'
#' @export
model_topographic <- function(noise_convention = c("split", "full"), params = NULL) {
  noise_convention <- match.arg(noise_convention)
  defaults <- list(gam = 0.18, nu = 0.02, nuV = 0.01,
                   s_psi = 0.03, s_V = 0.015)
  p <- .complete_scalar_params(params, defaults, "topographic",
    c("gam", "nu", "nuV", "s_psi", "s_V"))
  p$beta <- 5 * sqrt(p$gam)
  Kp <- rbind(c(0,2), c(1,2), c(0,1), c(1,1), c(2,1),
              c(1,0), c(2,0), c(1,-1), c(2,-1), c(1,-2))
  nK <- nrow(Kp); fullK <- rbind(Kp, -Kp)
  key <- function(k) paste(k[1], k[2], sep = "_")
  hpos <- c("1_1" = 0, "2_1" = -0.75, "1_0" = 0.15, "2_0" = 1.55,
            "1_-1" = 0.90, "2_-1" = 1.65) * p$gam
  hval <- function(k) { v <- hpos[key(k)]; if (is.na(v)) v <- hpos[key(-k)]
    if (is.na(v)) 0 else unname(v) }
  idx <- stats::setNames(seq_len(nK), apply(Kp, 1, key))
  phi_of <- function(z, k) {
    i <- idx[key(k)]; if (!is.na(i)) return(z[i])
    i <- idx[key(-k)]; if (!is.na(i)) return(Conj(z[i]))
    0 + 0i
  }
  in_K <- function(q) !is.na(idx[key(q)]) || !is.na(idx[key(-q)])
  triads <- lapply(seq_len(nK), function(i) {
    k <- Kp[i, ]; out <- list()
    for (r in seq_len(2 * nK)) {
      mv <- fullK[r, ]; q <- k - mv
      if (in_K(q)) {
        kperp <- c(-k[2], k[1]); mperp <- c(-mv[2], mv[1])
        out[[length(out) + 1]] <- list(m = mv, q = q,
          C = -sum((kperp - mperp) * mv), m2 = sum(mv^2), hm = hval(mv))
      }
    }
    out
  })
  k2 <- rowSums(Kp^2)
  drift_c <- function(t, V, z) {
    dz <- complex(nK)
    for (i in seq_len(nK)) {
      k <- Kp[i, ]; acc <- 0 + 0i
      for (tr in triads[[i]])
        acc <- acc + (tr$C / k2[i]) * phi_of(z, tr$q) *
                     (-tr$m2 * phi_of(z, tr$m) + tr$hm)
      dz[i] <- acc + 1i * k[1] * (p$beta / k2[i] - V) * z[i] +
               1i * k[1] * hval(k) / k2[i] * V - p$nu * z[i]
    }
    dV <- -p$nuV * V
    for (r in seq_len(2 * nK)) {
      kv <- fullK[r, ]
      dV <- dV + kv[1] * Im(hval(kv) * Conj(phi_of(z, kv)))
    }
    list(dz = dz, dV = dV)
  }
  sig_c <- if (noise_convention == "split") p$s_psi / sqrt(2) else p$s_psi
  m <- stochastic_model(
    f = function(t, x, y) { z <- complex(real = y[1:nK], imaginary = y[nK + 1:nK])
      drift_c(t, x[1], z)$dV },
    g = function(t, x, y) { z <- complex(real = y[1:nK], imaginary = y[nK + 1:nK])
      d <- drift_c(t, x[1], z)$dz; c(Re(d), Im(d)) },
    Sx = function(t, x) matrix(p$s_V, 1, 1),
    Sy = function(t, x, y) diag(sig_c, 2 * nK),
    k = 1, l = 2 * nK, name = "topographic (V | modes)")
  m$meta$params <- p
  mode_labels <- apply(Kp, 1, function(k) sprintf("%d_%d", k[1], k[2]))
  m$meta$vars <- list(observed = "V",
    hidden = c(paste0("Re_phi_", mode_labels),
               paste0("Im_phi_", mode_labels)))
  m$meta$provenance <- "moser2026extremes equations (4.5)-(4.13) and Appendix B"
  m$meta$source_status <- paste(
    "paper checked; no corresponding MATLAB supplied;",
    "complex-noise component convention remains an open discrepancy")
  m$meta$preset_caveat <- paste(
    "noise_convention='split' assigns sigma_psi/sqrt(2) to each real",
    "component; 'full' assigns sigma_psi to each. moser2026extremes does not disambiguate",
    "these complex-Wiener conventions.")
  m$meta$modes <- Kp
  m$meta$mode_families <- list(upper = c("1_1", "2_1"), zonal = c("1_0", "2_0"),
                               lower = c("1_-1", "2_-1"), merid = c("0_1", "0_2"))
  m$meta$noise_convention <- noise_convention
  m$meta$ic_default <- list(x0 = 0.5, y0 = rep(0.01, 2 * nK))
  m
}
