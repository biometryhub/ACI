# The Milstein correction for the wind-burst variable: an anomaly for review

**Date:** 13 August 2026
**Status:** Unresolved. Flagged for author review; aciR has neither replicated
nor silently corrected it.
**Reviewed at commit:** `7b25c66`
**Reference:** `matlab_reference/ENSO_model_cond_ACI_*.m`, upstream commit
`733c49fc5a36a5dad857608ef550a8f836466144`

---

## Summary

The reference implementation applies a Milstein correction to the intraseasonal
wind-burst variable `τ` using the derivative of `τ`'s diffusion with respect to
a **different variable** (`T_C`), multiplied against `τ`'s own squared Wiener
increment. Under the standard Milstein scheme that correction is zero, because
`σ_τ` does not depend on `τ`.

Either the reference term is an error, or it is a deliberate approximation
whose justification is not written down anywhere I can find. **I cannot
distinguish these from the available material**, so aciR implements the
standard scheme and this document records the discrepancy rather than
resolving it.

This matters for anyone reproducing the ENSO case study: the two choices give
different sample paths.

---

## What the reference does

From `ENSO_model_cond_ACI_u_h_W_tau_unobs.m`, the wind-burst update:

```matlab
% τ's multiplicative noise coefficient.
sigma_tau = 0.9 * (tanh(7.5*T_C(j-1)) + 1) * (1 + 0.3*cos((j-1)*dt*2*pi/6 + 2*pi/6));

% Derivative of τ's multiplicative noise coefficient.
der_sigma_tau = 0.9 * 7.5 * sech(7.5*T_C(j-1))^2 * (1 + 0.3*cos((j-1)*dt*2*pi/6 + 2*pi/6));

% Updating the intraseasonal wind burst, τ.
tau(j) = tau(j-1) - d_tau * tau(j-1) * dt ...
         + sigma_tau * sqrt(dt) * dW_tau ...
         + 0.5 * sigma_tau * der_sigma_tau * (dt * dW_tau^2 - dt);
```

Note that `der_sigma_tau` is `∂σ_τ/∂T_C` — the derivative with respect to the
central-Pacific sea-surface temperature, not with respect to `τ`.

## Why this looks wrong

For a scalar stochastic differential equation

```
dX = a(X) dt + b(X) dW
```

the Milstein scheme is

```
X_{n+1} = X_n + a Δt + b ΔW + ½ b (∂b/∂X) (ΔW² − Δt)
```

The correction term carries the derivative of the diffusion **with respect to
the variable being integrated**, and it exists because `b` changes over the
step as `X` itself moves.

Here `X = τ` and `b = σ_τ(t, T_C)`. Since `σ_τ` does not depend on `τ`:

```
∂σ_τ/∂τ = 0   ⟹   the Milstein correction for τ is zero
```

and the scheme reduces to Euler–Maruyama, exactly as the reference itself
notes for the four interannual variables whose "noise feedbacks are constant or
state-independent".

The reference instead uses `∂σ_τ/∂T_C` against `dW_τ²`. That mixes the
derivative with respect to one variable with the squared increment of a
different, **independent** Wiener process. A correction accounting for `T_C`'s
motion over the step would involve the cross term `∂σ_τ/∂T_C · σ_C` and the
product of increments `ΔW_τ ΔW_C`, not `ΔW_τ²`, and it would not be the
Milstein scheme for a scalar equation.

## Contrast: the Walker circulation is correct

The same file's treatment of `I` is textbook, which makes the `τ` term stand
out rather than look like a house convention:

```matlab
sigma_I = sqrt(lambda * (4-I(j-1)) * I(j-1));
der_sigma_I = lambda * (2-I(j-1)) / sigma_I;
I(j) = I(j-1) - lambda * (I(j-1) - m) * dt + sigma_I * sqrt(dt) * dW_I ...
       + 0.5 * sigma_I * der_sigma_I * (dt * dW_I^2 - dt);
```

Here `σ_I` genuinely depends on `I`, the variable being integrated,
`der_sigma_I` is `∂σ_I/∂I`, and the correction is the standard one. The
comment even explains why it is needed and what equilibrium distribution it
preserves.

## What I checked before flagging this

- **The Nature Communications paper**: 9 pages, main article only. No mention.
- **The full preprint with Supplementary Information**
  (arXiv 2505.14825, 47 pages, SI on 18–47): the word "Milstein" **does not
  appear anywhere in the document.** The integration scheme is an
  implementation detail of the MATLAB, not a documented part of the method.
- **The online-smoother paper** (arXiv 2411.05870): concerns the estimator, not
  the forward simulation.

So the paper cannot adjudicate this, and I have no basis to prefer one reading
over the other beyond the mathematics.

## What aciR does

`aci_simulate(scheme = "milstein")` implements the **standard** correction:
the derivative is taken with respect to the integrated variable, and the caller
supplies it explicitly (`d_sigma_x`) rather than having it estimated. It is
graded by strong convergence order — measured 1.01, 1.14, 0.96 against the
theoretical 1.0, on geometric Brownian motion driven by the same Wiener path
the simulator integrated, with Euler–Maruyama measuring 0.66, 0.47, 0.55
against its theoretical 0.5.

aciR therefore does **not** reproduce the reference's `τ` update. Any
comparison of simulated ENSO wind-burst paths between the two implementations
will differ, and the difference is this term.

## The three possibilities, and what would settle each

| Possibility | What would settle it |
|---|---|
| The reference term is an error | Confirmation from the authors, or a derivation showing the standard term is intended |
| It is a deliberate approximation | A stated justification — e.g. treating `T_C`'s motion as locally driven by the same increment |
| I have misread the intended SDE | A statement of the `τ` equation in which `σ_τ` depends on `τ` |

## Why it is left open rather than decided

Replicating a term I believe to be wrong would put a defect into aciR under the
banner of fidelity. Silently "correcting" the authors' published code would
mean aciR diverges from the reference at a point no reader is told about. Both
are worse than saying plainly that the two differ here and why.

The decision to raise this with the authors is **Max's alone**; nothing in this
document has been sent anywhere.

## Practical consequence

Nothing currently shipped in aciR depends on this. The Milstein scheme is a
simulation capability, `aci_enso_model()` supplies components from observed
paths rather than simulating them, and no shipped fixture exercises the `τ`
update. The anomaly becomes load-bearing only if and when the ENSO case study
is reproduced end to end, which is where it should be settled.
