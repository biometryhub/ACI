## acir reserve file
## Origin: aci/R/assimilation.R:709-712
## Source package: aci 0.0.30, git tree 97f6b124
## Category: fbcir
## Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
## Reason: Unreachable once mode='one_lag' is gone.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The is.null(x$rows) arm of as.data.frame.lag_table, unreachable once mode='one_lag' is gone.

  if (is.null(x$rows)) {
    ol <- lt_onelag(x)
    return(data.frame(j = seq_along(ol), t = x$t[seq_along(ol)], onelag = ol))
  }
