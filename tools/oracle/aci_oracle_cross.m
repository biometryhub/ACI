% aci_oracle_cross.m, cross-noise oracle harness for `oracle.aci_matlab_reference`.
%
% WHY THIS EXISTS. The dyad oracle (aci_oracle_dyad.m) grades aciR's core against
% the authors' reference, but that reference, like EVERY scalar model upstream
% (dyad_interaction_model.m, noisy_predator_prey_model.m), sets Sx_2 = 0 and
% Sy_1 = 0, hence S_yoS_x = 0. The noise-cross-covariance terms of the filter
% (aux = S_yoS_x + R*L_x) and smoother (A_j, and the S_yoS_x transport term)
% therefore EXECUTE on every dyad run but are annihilated: the fixture never
% grades them. aciR exposes those paths publicly via aci_cgns_model(S_yoS_x = ..),
% so they need an oracle of their own.
%
% WHAT THIS GRADES, AND WHAT IT DOES NOT. This harness is a second, independent
% implementation of the published CGNS equations in a different language and
% runtime; agreement with it refutes an R-side transcription error on the full
% transient. It is NOT an authors'-reference grounding of the cross-noise terms,
% because upstream has no scalar model that exercises them. The primary grounding
% for those terms is analytic, the algebraic Kalman-Bucy fixed points, solved in
% closed form in the package's tests/testthat/test-identities.R, independent of
% BOTH implementations. The two are complementary: analytic covers the stationary
% regime exactly, this covers the transient.
%
% The system is the reference dyad with correlated noise switched on, so the only
% difference from the graded flagship is the term under test.
% Run: matlab -nodisplay -nosplash -batch "run('aci_oracle_cross.m')"
% Generated on MATLAB R2025b Update 1 (25.2.0.3042426).
%
% Portions derived from dyad_interaction_model.m (github.com/marandmath/ACI_code),
% Copyright (c) 2025 Marios Andreou, released under the MIT License. That
% copyright notice and the permission notice in matlab_reference/LICENSE apply
% to the derived portions.

rng(4242)                                  % distinct from the dyad harness's 333
N = 30000; dt = 0.001;

% model + parameters (the reference dyad; noise decomposition changed) ---------
x = zeros(1, N+1); y = zeros(1, N+1);
d_x = 0.5; d_y = 0.5; gamma = 2; F_x = 0.5; F_y = 1;
x(1) = F_x/d_x; y(1) = F_y/d_y;
L_x = zeros(1, N+1); f_x = zeros(1, N+1);
% Correlated noise: each process is driven by BOTH Wiener increments.
%   S_xoS_x = Sx_1^2 + Sx_2^2 = 0.45
%   S_yoS_y = Sy_1^2 + Sy_2^2 = 0.89
%   S_yoS_x = Sy_1*Sx_1 + Sy_2*Sx_2 = 0.54   (non-zero, the term under test)
%   det = S_xoS_x*S_yoS_y - S_yoS_x^2 = 0.1089 > 0, so the joint covariance is
%   positive definite and the system is admissible.
Sx_1 = 0.6; Sx_2 = 0.3;
L_y = -d_y; f_y = zeros(1, N+1); Sy_1 = 0.5; Sy_2 = 0.8;
L_x(1) = gamma*x(1); f_x(1) = F_x - d_x*x(1); f_y(1) = F_y - gamma*x(1)^2;

% true signals -----------------------------------------------------------------
for j = 2:N+1
    dW_x = randn; dW_y = randn;
    x(j) = x(j-1) + (L_x(j-1)*y(j-1) + f_x(j-1))*dt + Sx_1*sqrt(dt)*dW_x + Sx_2*sqrt(dt)*dW_y;
    y(j) = y(j-1) + (L_y*y(j-1) + f_y(j-1))*dt + Sy_1*sqrt(dt)*dW_x + Sy_2*sqrt(dt)*dW_y;
    L_x(j) = gamma*x(j); f_x(j) = F_x - d_x*x(j); f_y(j) = F_y - gamma*x(j)^2;
end

% CGNS filter (forward) --------------------------------------------------------
S_xoS_x = Sx_1^2 + Sx_2^2; S_xoS_x_inv = 1/S_xoS_x;
S_yoS_y = Sy_1^2 + Sy_2^2; S_yoS_x = Sy_1*Sx_1 + Sy_2*Sx_2; S_xoS_y = S_yoS_x.';
filter_mean = zeros(1, N+1); filter_mean(1) = y(1);
filter_cov  = zeros(1, N+1); filter_cov(1) = 0.1;
mu0 = filter_mean(1); R0 = filter_cov(1);
for j = 2:N+1
    dx = x(j) - x(j-1); aux = S_yoS_x + filter_cov(j-1)*L_x(j-1);
    mu = mu0 + (L_y*mu0 + f_y(j-1))*dt + aux*S_xoS_x_inv*(dx - (L_x(j-1)*mu0 + f_x(j-1))*dt);
    R  = R0  + (L_y*R0 + R0*L_y + S_yoS_y - aux*S_xoS_x_inv*aux)*dt;
    filter_mean(j) = mu; filter_cov(j) = R; mu0 = mu; R0 = R;
end

% CGNS smoother (backward) -----------------------------------------------------
smoother_mean = zeros(1, N+1); smoother_cov = zeros(1, N+1);
smoother_mean(N+1) = filter_mean(N+1); smoother_cov(N+1) = filter_cov(N+1);
muT = smoother_mean(N+1); RT = smoother_cov(N+1);
for j = N:-1:1
    dx = x(j+1) - x(j);
    A_j = L_y - S_yoS_x*S_xoS_x_inv*L_x(j);
    B_j = S_yoS_y - S_yoS_x*S_xoS_x_inv*S_xoS_y;
    mu = muT - (L_y*muT + f_y(j) - B_j/filter_cov(j)*(filter_mean(j) - muT))*dt ...
             + S_yoS_x*S_xoS_x_inv*(-dx + (L_x(j)*muT + f_x(j))*dt);
    R  = RT - ((A_j + B_j/filter_cov(j))*RT + RT*(A_j + B_j/filter_cov(j)) - B_j)*dt;
    smoother_mean(j) = mu; smoother_cov(j) = R; muT = mu; RT = R;
end

% ACI metric (relative entropy: smoother || filter) ----------------------------
signal = 0.5*(smoother_mean - filter_mean).^2 ./ filter_cov;
cov_ratio = smoother_cov ./ filter_cov;
dispersion = 0.5*(-log(cov_ratio) + cov_ratio - 1);
ACI_metric = signal + dispersion;

% guard: the fixture is worthless if the run left the admissible domain --------
assert(all(isfinite(filter_cov)) && all(filter_cov > 0), 'filter covariance left its domain');
assert(all(isfinite(smoother_cov)) && all(smoother_cov > 0), 'smoother covariance left its domain');
assert(S_yoS_x ~= 0, 'this harness must exercise a non-zero cross-covariance');

% write fixtures ---------------------------------------------------------------
writematrix([(0:N)'*dt, x'], 'cross_signal_x.csv');
idx = 1:100:N+1; t = (idx-1)*dt;
Tbl = table(t(:), x(idx)', y(idx)', filter_mean(idx)', filter_cov(idx)', ...
            smoother_mean(idx)', smoother_cov(idx)', ACI_metric(idx)', ...
    'VariableNames', {'t','x','y','filter_mean','filter_cov', ...
                      'smoother_mean','smoother_cov','ACI_metric'});
writetable(Tbl, 'cross_reference.csv');

% deterministic checksums (the oracle's staleness/identity anchor)
fprintf('OK rows=%d N=%d\n', numel(idx), N);
fprintf('GRAMMIANS S_xoS_x=%.10f S_yoS_y=%.10f S_yoS_x=%.10f det=%.10f\n', ...
        S_xoS_x, S_yoS_y, S_yoS_x, S_xoS_x*S_yoS_y - S_yoS_x^2);
fprintf('CHECKSUM sum(x)=%.10f sum(filter_mean)=%.10f sum(smoother_mean)=%.10f sum(ACI_metric)=%.10f\n', ...
        sum(x), sum(filter_mean), sum(smoother_mean), sum(ACI_metric));
fprintf('SPOT x(30001)=%.10f filter_mean(30001)=%.10f ACI_metric(15001)=%.10f\n', ...
        x(N+1), filter_mean(N+1), ACI_metric(15001));
