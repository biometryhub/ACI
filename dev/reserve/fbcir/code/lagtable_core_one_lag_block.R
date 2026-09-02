## acir reserve file
## Origin: aci/R/compiled_lag.R:126-156
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: The mode=='one_lag' early return of the compiled lag core.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The mode=='one_lag' early-return block of .lagtable_core_compiled; self-contained, ends with its own return().

  if (mode == "one_lag") {
    coN1 <- .compiled_co(bundle, N)
    auxN1 <- .thmD1_aux_compiled(bundle, N, CVf[, , N], co = coN1)
    ol <- .onelag_stats(
      coN1, auxN1, MUf[N, ], CVf[, , N], MUf[N1, ], CVf[, , N1],
      bundle$x[N1, ] - bundle$x[N, ], dt, l
    )
    dmu <- ol$mu - MUf[N, ]
    dR <- ol$R - CVf[, , N]
    P <- rep(0, N1)
    D <- diag(l)
    for (jj in N:1) {
      if (jj < N) {
        co <- .compiled_co(bundle, jj)
        E <- .thmD1_aux_compiled(bundle, jj, CVf[, , jj], co = co)$E
        D <- E %*% D
      }
      RsjN <- smoo$cov[, , jj]
      A <- sym(D %*% dR %*% t(D))
      Rlag <- spd_floor(RsjN - A)
      v <- drop(D %*% dmu)
      sig <- 0.5 * sum(v * chol_solve(Rlag, v, "Rlag"))
      trA <- sum(diag(chol_solve(Rlag, A, "Rlag")))
      ld <- logdet_chol(RsjN) - logdet_chol(Rlag)
      P[jj] <- max(sig + 0.5 * (trA - ld), 0)
    }
    return(list(
      diag = diagv, rows = NULL, L = NULL, tailbnd = NULL,
      onelag = P, stop_index = NA_integer_
    ))
  }
