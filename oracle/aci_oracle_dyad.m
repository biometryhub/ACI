% aci_oracle_dyad.m — foundry oracle harness for `oracle.aci_matlab_reference`.
%
% Reproduces the DETERMINISTIC ACI core (CGNS filter -> smoother -> ACI metric) of
% dyad_interaction_model.m (github.com/marandmath/ACI_code, MIT, Andreou/Chen/Bollt)
% with the reference's exact rng(333) signals, and writes two portable CSV fixtures:
%   - dyad_signal_x.csv    : the full observable signal x (the INPUT aciR runs on)
%   - dyad_reference.csv   : the expected filter/smoother/ACI outputs (subsampled)
% aciR is graded (IOP) by running ITS OWN filter/smoother/metric on dyad_signal_x
% and matching dyad_reference at the sampled indices to numerical tolerance.
%
% Core only (O(N), seconds); the O(N^2) causal-influence-range is a later fixture.
% Run: matlab -batch "run('aci_oracle_dyad.m')"   (tested target MATLAB R2024b/R2025b)
%
% Portions derived from dyad_interaction_model.m, Copyright (c) 2025 Marios
% Andreou, released under the MIT License. That copyright notice and the
% permission notice in matlab_reference/LICENSE apply to the derived portions.

rng(333)                                   % the reference's seed — authentic signals
N = 30000; dt = 0.001;

% -- model + parameters (verbatim from dyad_interaction_model.m) ---------------
x = zeros(1, N+1); y = zeros(1, N+1);
d_x = 0.5; d_y = 0.5; gamma = 2; F_x = 0.5; F_y = 1; sigma_x = 0.5; sigma_y = 1;
x(1) = F_x/d_x; y(1) = F_y/d_y;
L_x = zeros(1, N+1); f_x = zeros(1, N+1);
Sx_1 = sigma_x; Sx_2 = 0;
L_y = -d_y; f_y = zeros(1, N+1); Sy_1 = 0; Sy_2 = sigma_y;
L_x(1) = gamma*x(1); f_x(1) = F_x - d_x*x(1); f_y(1) = F_y - gamma*x(1)^2;

% -- true signals -------------------------------------------------------------
for j = 2:N+1
    dW_x = randn; dW_y = randn;
    x(j) = x(j-1) + (L_x(j-1)*y(j-1) + f_x(j-1))*dt + Sx_1*sqrt(dt)*dW_x + Sx_2*sqrt(dt)*dW_y;
    y(j) = y(j-1) + (L_y*y(j-1) + f_y(j-1))*dt + Sy_1*sqrt(dt)*dW_x + Sy_2*sqrt(dt)*dW_y;
    L_x(j) = gamma*x(j); f_x(j) = F_x - d_x*x(j); f_y(j) = F_y - gamma*x(j)^2;
end

% -- CGNS filter (forward) ----------------------------------------------------
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

% -- CGNS smoother (backward) -------------------------------------------------
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

% -- ACI metric (relative entropy: smoother || filter) ------------------------
signal = 0.5*(smoother_mean - filter_mean).^2 ./ filter_cov;
cov_ratio = smoother_cov ./ filter_cov;
dispersion = 0.5*(-log(cov_ratio) + cov_ratio - 1);
ACI_metric = signal + dispersion;

% -- write fixtures -----------------------------------------------------------
% full input signal x (aciR runs its own filter/smoother on THIS)
writematrix([(0:N)'*dt, x'], 'dyad_signal_x.csv');
% expected outputs, subsampled every 100th point (301 rows) for a compact fixture
idx = 1:100:N+1; t = (idx-1)*dt;
Tbl = table(t(:), x(idx)', y(idx)', filter_mean(idx)', filter_cov(idx)', ...
            smoother_mean(idx)', smoother_cov(idx)', ACI_metric(idx)', ...
    'VariableNames', {'t','x','y','filter_mean','filter_cov', ...
                      'smoother_mean','smoother_cov','ACI_metric'});
writetable(Tbl, 'dyad_reference.csv');

% deterministic checksums (the oracle's staleness/identity anchor)
fprintf('OK rows=%d N=%d\n', numel(idx), N);
fprintf('CHECKSUM sum(x)=%.10f sum(filter_mean)=%.10f sum(smoother_mean)=%.10f sum(ACI_metric)=%.10f\n', ...
        sum(x), sum(filter_mean), sum(smoother_mean), sum(ACI_metric));
fprintf('SPOT x(30001)=%.10f filter_mean(30001)=%.10f ACI_metric(15001)=%.10f\n', ...
        x(N+1), filter_mean(N+1), ACI_metric(15001));
