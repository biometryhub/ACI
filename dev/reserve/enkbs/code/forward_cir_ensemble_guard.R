## acir reserve file
## Origin: aci/R/causal_metrics.R:766-769
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: Ensemble guard at the head of forward_cir.aci_result; unreachable without the ensemble engine.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The ensemble guard abort at the head of forward_cir.aci_result (cut point 3); unreachable once the ensemble engine is out.

  if (is.null(x$table) && identical(x$meta$engine, "ensemble"))
    aci_abort("aci_error_not_implemented", paste(
      "The forward ensemble CIR needs the lagged EnKBS family; recompute",
      "aci(..., engine = 'ensemble', keep = 'table')."))
