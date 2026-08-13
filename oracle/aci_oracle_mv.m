% ACI_ORACLE_MV  Independent-oracle fixture for the VECTOR-valued CGNS core,
% with a non-zero matrix noise cross-covariance.
%
% THIS HARNESS HAS NO UPSTREAM COUNTERPART, AND THE GAP IS LARGER HERE THAN
% ANYWHERE ELSE IN THE PACKAGE.
%
% Every scalar model in the reference implementation sets the noise
% cross-covariance to zero, and every ENSO script says so in as many words --
% "NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT DEFINE
% THE ASSOCIATED MATRIX-VALUED FUNCTIONALS IN THIS SCRIPT FOR SIMPLICITY",
% thirty-five times across the five files. So the matrix-valued cross-noise
% terms that aciR exposes publicly are graded by NOTHING upstream.
%
% This file is therefore a SECOND INDEPENDENT IMPLEMENTATION of the published
% matrix CGNS equations, in a different language and runtime. Agreement refutes
% an R-side transcription error. It is NOT an authors'-reference grounding, and
% the manifest must not claim it is. The primary grounding for these terms is
% analytic: the block-diagonal algebraic Riccati fixed points solved in closed
% form in tests/testthat/test-identities-mv.R, which depend on neither
% implementation.
%
% Run from the oracle/ directory:
%   matlab -nodisplay -nosplash -batch "run('aci_oracle_mv.m')"

clearvars; close all;
rng(2718)

N = 4000;
dt = 0.002;
n_x = 2;
n_y = 2;

% ---- A two-by-two system with genuinely coupled, correlated noise -----------
% The noise decomposition is written as a single joint factor S so that the
% joint covariance is positive semidefinite by construction rather than by
% assertion, and the cross-block is non-zero by construction too.
S = [ 0.60  0.20  0.15  0.05;
      0.10  0.50  0.05  0.10;
      0.25  0.10  0.70  0.15;
      0.05  0.30  0.10  0.55];
S_x = S(1:n_x, :);
S_y = S(n_x+1:end, :);
S_xoS_x = S_x * S_x';
S_yoS_y = S_y * S_y';
S_yoS_x = S_y * S_x';
S_xoS_x_inv = inv(S_xoS_x);

fprintf('CROSSNOISE_FROBENIUS %.6e\n', norm(S_yoS_x, 'fro'));

% ---- Time-varying coefficients ----------------------------------------------
L_x = zeros(n_x, n_y, N+1);
L_y = zeros(n_y, n_y, N+1);
f_x = zeros(n_x, N+1);
f_y = zeros(n_y, N+1);

x = zeros(n_x, N+1);
y = zeros(n_y, N+1);
x(:, 1) = [0.5; -0.3];
y(:, 1) = [1.0;  0.4];

    function [Lx, Ly, fx, fy] = coeffs(xv, t)
        Lx = [ 0.8 + 0.3*xv(1),      0.2*sin(t);
               0.1*xv(2),            0.6 - 0.2*xv(1) ];
        Ly = [-1.2 + 0.1*xv(1),      0.3;
               0.2,                 -0.9 - 0.1*xv(2) ];
        fx = [ 0.4 - 0.5*xv(1); -0.3*xv(2) + 0.2 ];
        fy = [ 0.5 - 0.2*xv(1)^2;  0.1 - 0.15*xv(2) ];
    end

[Lx0, Ly0, fx0, fy0] = coeffs(x(:,1), 0);
L_x(:,:,1) = Lx0; L_y(:,:,1) = Ly0; f_x(:,1) = fx0; f_y(:,1) = fy0;

for j = 2:N+1
    dW = randn(4, 1) * sqrt(dt);
    x(:, j) = x(:, j-1) + (L_x(:,:,j-1) * y(:, j-1) + f_x(:, j-1)) * dt ...
              + S_x * dW;
    y(:, j) = y(:, j-1) + (L_y(:,:,j-1) * y(:, j-1) + f_y(:, j-1)) * dt ...
              + S_y * dW;
    [Lxj, Lyj, fxj, fyj] = coeffs(x(:, j), (j-1)*dt);
    L_x(:,:,j) = Lxj; L_y(:,:,j) = Lyj; f_x(:,j) = fxj; f_y(:,j) = fyj;
end

% ---- Filter -----------------------------------------------------------------
filter_mean = zeros(n_y, N+1);
filter_cov = zeros(n_y, n_y, N+1);
filter_mean(:, 1) = [0.8; 0.2];
filter_cov(:, :, 1) = 0.2 * eye(n_y);
mu = filter_mean(:, 1);
R = filter_cov(:, :, 1);

for j = 2:N+1
    dx = x(:, j) - x(:, j-1);
    aux = S_yoS_x + R * L_x(:,:,j-1)';
    mu = mu + (L_y(:,:,j-1) * mu + f_y(:, j-1)) * dt ...
            + aux * S_xoS_x_inv * (dx - (L_x(:,:,j-1) * mu + f_x(:, j-1)) * dt);
    R = R + (L_y(:,:,j-1) * R + R * L_y(:,:,j-1)' + S_yoS_y ...
             - aux * S_xoS_x_inv * aux') * dt;
    R = (R + R') / 2;
    filter_mean(:, j) = mu;
    filter_cov(:, :, j) = R;
end

% ---- Smoother ---------------------------------------------------------------
smoother_mean = zeros(n_y, N+1);
smoother_cov = zeros(n_y, n_y, N+1);
smoother_mean(:, N+1) = filter_mean(:, N+1);
smoother_cov(:, :, N+1) = filter_cov(:, :, N+1);
muT = smoother_mean(:, N+1);
RT = smoother_cov(:, :, N+1);
S_xoS_y = S_yoS_x';

for j = N:-1:1
    dx = x(:, j+1) - x(:, j);
    A_j = L_y(:,:,j) - S_yoS_x * S_xoS_x_inv * L_x(:,:,j);
    B_j = S_yoS_y - S_yoS_x * S_xoS_x_inv * S_xoS_y;
    Rf = filter_cov(:, :, j);
    BRf = B_j / Rf;
    mu = muT - (L_y(:,:,j) * muT + f_y(:, j) ...
                - BRf * (filter_mean(:, j) - muT)) * dt ...
             + S_yoS_x * S_xoS_x_inv * (-dx + (L_x(:,:,j) * muT + f_x(:, j)) * dt);
    M = A_j + BRf;
    R = RT - (M * RT + RT * M' - B_j) * dt;
    R = (R + R') / 2;
    smoother_mean(:, j) = mu;
    smoother_cov(:, :, j) = R;
    muT = mu; RT = R;
end

% ---- Multivariate relative entropy -----------------------------------------
aci = zeros(1, N+1);
for j = 1:N+1
    Rf = filter_cov(:, :, j);
    Rs = smoother_cov(:, :, j);
    d = smoother_mean(:, j) - filter_mean(:, j);
    aci(j) = 0.5 * (d' / Rf * d + trace(Rf \ Rs) - n_y ...
                    - log(det(Rs) / det(Rf)));
end

% ---- Write fixtures ---------------------------------------------------------
idx = 1:20:(N+1);
signal = array2table([(0:N)' * dt, x'], 'VariableNames', {'t', 'x1', 'x2'});
writetable(signal, 'mv_signal.csv');

ref = array2table([idx', ...
    filter_mean(1, idx)', filter_mean(2, idx)', ...
    squeeze(filter_cov(1,1,idx)), squeeze(filter_cov(1,2,idx)), squeeze(filter_cov(2,2,idx)), ...
    smoother_mean(1, idx)', smoother_mean(2, idx)', ...
    squeeze(smoother_cov(1,1,idx)), squeeze(smoother_cov(1,2,idx)), squeeze(smoother_cov(2,2,idx)), ...
    aci(idx)'], ...
    'VariableNames', {'index','fm1','fm2','fc11','fc12','fc22', ...
                      'sm1','sm2','sc11','sc12','sc22','ACI_metric'});
writetable(ref, 'mv_reference.csv');

fprintf('MV_ORACLE_DONE rows=%d metric_range=[%.6f %.6f]\n', ...
        numel(idx), min(aci), max(aci));
