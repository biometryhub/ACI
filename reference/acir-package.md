# acir: Closed-Form Assimilative Causal Inference for Conditional Gaussian Nonlinear Systems

Causal inference as a Bayesian inverse problem via data assimilation,
restricted to the closed-form conditional-Gaussian (CGNS) engine.
Implements the assimilative causal inference (ACI) metric, closed-form
CGNS filtering and smoothing, the Theorem 3 lagged-divergence table with
its online smoother, conditional ACI masking, the forward causal
influence range, and the benchmark models of the reference MATLAB
codebase. Based on Andreou, Chen & Bollt (2026)
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)
, with the causal influence range of Andreou & Chen (2026)
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)
and the adaptive online smoother of Andreou, Chen & Li (2026)
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)
.

## See also

Useful links:

- <https://github.com/biometryhub/ACI>

- <https://biometryhub.github.io/ACI/>

- Report bugs at <https://github.com/biometryhub/ACI/issues>

## Author

**Maintainer**: Aidan Moller <aidan.moller@adelaide.edu.au> \[copyright
holder\]

Authors:

- Max Moldovan ([ORCID](https://orcid.org/0000-0001-9680-8474))
  \[copyright holder\]

Other contributors:

- Marios Andreou (Author and copyright holder of the MIT-licensed
  reference MATLAB implementation 'ACI_code', against which this package
  is verified and from which parts of its numerical core derive; see
  inst/COPYRIGHTS) \[copyright holder\]
