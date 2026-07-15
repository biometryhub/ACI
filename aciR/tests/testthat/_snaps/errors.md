# the model constructor's error messages are stable

    Code
      aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1, S_yoS_y = -1)
    Condition
      Error:
      ! `S_yoS_y` must be non-negative; it is the latent-noise covariance.

---

    Code
      aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1, S_yoS_y = 1,
        S_yoS_x = 1.5)
    Condition
      Error:
      ! The joint noise covariance is not positive semidefinite: `S_xoS_x` * `S_yoS_y` - `S_yoS_x`^2 is -1.25. The noise cross-covariance cannot exceed sqrt(`S_xoS_x` * `S_yoS_y`) = 1 in magnitude.

---

    Code
      aci_cgns_model(L_x = 1, f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 0, S_yoS_y = 1)
    Condition
      Error:
      ! `S_xoS_x` must be positive; it is the observation-noise covariance.

---

    Code
      aci_cgns_model(L_x = "nope", f_x = 0, L_y = -0.5, f_y = 0, S_xoS_x = 1,
        S_yoS_y = 1)
    Condition
      Error:
      ! `L_x` must be a vectorised function of the observed signal or a finite numeric scalar; it is character of length 1.

# the observed-signal error messages are stable

    Code
      aci(c(1, NA, 3), model)
    Condition
      Error:
      ! `x` must be complete and finite; it holds NA at index 2. The filter has no missing-observation contract: interpolate or subset to a complete, regularly sampled span before calling.

---

    Code
      aci(c(1, Inf, 3), model)
    Condition
      Error:
      ! `x` must be complete and finite; it holds an infinite value at index 2. The filter has no missing-observation contract: interpolate or subset to a complete, regularly sampled span before calling.

---

    Code
      aci(1, model)
    Condition
      Error:
      ! `x` must be a numeric vector of at least two observations.

---

    Code
      aci(matrix(1:10, nrow = 2L), model)
    Condition
      Error:
      ! `x` must be a plain numeric vector; it is an object of class matrix with dimension 2 x 5.

# the coefficient-contract error messages are stable

    Code
      aci(c(1, 2, 3, 4), scalar_coef)
    Condition
      Error:
      ! `L_x` must return one value per observation: expected length 4, got 1. Coefficient functions are not recycled; supply a vectorised function of the observed signal, or a numeric scalar for a coefficient that is constant in time.

# the components-schema error messages are stable

    Code
      aci_filter(x, comp[-1L], 0.001, 2, 0.1)
    Condition
      Error:
      ! `comp` is missing the component `L_x`. See `?aci_components` for the expected shape.

---

    Code
      aci_filter(x, asymmetric, 0.001, 2, 0.1)
    Condition
      Error:
      ! `comp$S_xoS_y` must equal `comp$S_yoS_x`: the noise cross-covariance of a scalar system is symmetric, but they are 0.2 and 0.1.

# the posterior-contract error messages are stable

    Code
      aci_metric(filt, list(mean = c(1, 2), cov = c(1, 1)))
    Condition
      Error:
      ! `smooth` must cover 3 time steps; it covers 2.

---

    Code
      aci_metric(filt, list(mean = c(1, 2, 3), cov = c(1, 0, 1)))
    Condition
      Error:
      ! `smooth$cov` must be finite and strictly positive; it holds 0 at index 2. A non-positive posterior covariance has no relative-entropy interpretation.

# the runtime covariance-guard message is stable

    Code
      aci(x, model, dt = 1)
    Condition
      Error:
      ! The filter covariance must stay finite and strictly positive; it became -10.15666 at index 3 (time 2). Reduce `dt`, or check the model's noise covariance: an explicit Euler step too large for the system can drive the covariance out of its domain even when the model is admissible.

# the time-grid error messages are stable

    Code
      aci(x, model, time = c(0, 1, 2, 4, 5))
    Condition
      Error:
      ! `time` must be equally spaced: the recursions integrate a fixed step, and the spacing ranges from 1 to 2. Resample onto a regular grid before calling.

---

    Code
      aci(x, model, time = c(0, 1, 2, 2, 3))
    Condition
      Error:
      ! `time` must be strictly increasing; it is not at index 4.

