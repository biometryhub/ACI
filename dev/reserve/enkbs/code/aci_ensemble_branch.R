## acir reserve file
## Origin: aci/R/causal_metrics.R:383-415
## Source package: aci 0.0.30, git tree 97f6b124
## Category: enkbs
## Intended release: 0.2.0 or 0.3.0 (EnKBS_code family; order TBD)
## Reason: The engine=='ensemble' arm of aci(); the only engine-side call into the excluded set.
## Verbatim excision from the acir 0.1.0 mainline from the aci 0.0.30 sources; not modified.
## Note: The engine=='ensemble' arm of aci() plus the engine-conditional KL path and lag-table calls it guarded (cut points 1-2). Retained here verbatim including the two closed-form lines interleaved with it.

  } else {
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
        m != floor(m) || m < 2L)
      aci_abort("aci_error_ensemble_rank", "m must be an integer of at least 2.")
    if (m <= model$l)
      aci_abort("aci_error_ensemble_rank", paste(
        "Full-dimensional Gaussian ACI from ensemble moments requires m > l;",
        "localized EnKBF/EnKBS state estimation may use m <= l, but its raw",
        "sample covariance is singular and needs a separately specified",
        "projection or covariance regularizer for KL."))
    fr <- enkbf(model, obs, m = m, seed = seed, nontarget = nontarget,
                localization = localization, inflation = inflation,
                ic_sampler = ic_sampler, ...)
    sm <- enkbs(fr$model, fr$path, fr$noise, localization = localization)
    filt <- as_gaussian(fr$path); smoo <- as_gaussian(sm)
    aci_warn("aci_warn_ensemble_kl", sprintf(
      "ACI from ensemble moments (m = %d): jiang2026enkbs s3.2 finds m ~ 50 faithful; m = 10 preserves timing/sign but distorts magnitudes.", m))
  }
  klp <- if (engine == "cgns")
    .gaussian_kl_path_compiled(
      bundle, smoo, filt, decompose = decompose, validate = FALSE
    ) else gaussian_kl_path(smoo, filt, decompose = decompose)
  tab <- NULL
  if (keep == "table") {
    tab <- if (engine == "cgns")
      .lag_table_compiled(
        bundle, mode = "forward", filter = filt, init = filt$meta$init,
        validate = FALSE
      ) else
      .ensemble_lag_table_from_run(model, obs, fr, sm,
                                   nontarget = nontarget,
                                   localization = localization)
  }
