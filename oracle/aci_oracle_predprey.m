% ACI_ORACLE_PREDPREY  Independent-oracle fixture for the predator-prey model.
%
% Transcribes the noisy predator-prey sections of the reference implementation
% (noisy_predator_prey_model.m, upstream commit
% 733c49fc5a36a5dad857608ef550a8f836466144) and writes the reference values
% aciR is graded against.
%
% This is the model that exercises a TIME-VARYING latent self-drift: in both
% causal directions the latent variable's own damping is set by the observed
% state. The dyad fixture cannot grade that path, because the dyad's latent
% self-drift is the constant -d_y.
%
% The reference studies the system in two directions and expects one section
% to be run at a time, each after rng(42). Both directions integrate the same
% Lotka-Volterra pair with the same noise applied the same way -- only the
% bookkeeping names swap -- so a single simulation serves both, and this
% harness verifies that equivalence rather than assuming it.
%
% Run from the oracle/ directory:
%   matlab -nodisplay -nosplash -batch "run('aci_oracle_predprey.m')"

clearvars; close all;

% ---- Parameters of the reference predator-prey model ------------------------
N = 12000;
dt = 0.005;
alpha = 0.4;       % predator's natural death rate
beta = 0.1;        % effect of prey availability on predator growth
gamma = 1.1;       % prey's natural growth rate
delta = 0.4;       % effect of predator presence on prey decline
sigma_x = 0.3;     % predator noise
sigma_y = 0.3;     % prey noise

% ---- Simulate the true signals (direction x -> y formulation) ---------------
rng(42)
x = zeros(1, N+1);      % predator
y = zeros(1, N+1);      % prey
x(1) = 4;
y(1) = 4;
L_x = zeros(1, N+1);    % unobservable coefficient: feedback of x in x
f_x = 0;
Sx_1 = sigma_x; Sx_2 = 0;
L_y = zeros(1, N+1);    % observable coefficient: feedback of x in y
f_y = zeros(1, N+1);
Sy_1 = 0; Sy_2 = sigma_y;
L_y(1) = - delta * y(1);
f_y(1) = gamma * y(1);
L_x(1) = beta * y(1) - alpha;

for j = 2:N+1
    dW_x = randn;
    dW_y = randn;
    x(j) = x(j-1) + (L_x(j-1) * x(j-1) + f_x) * dt ...
                  + Sx_1 * sqrt(dt) * dW_x + Sx_2 * sqrt(dt) * dW_y;
    y(j) = y(j-1) + (L_y(j-1) * x(j-1) + f_y(j-1)) * dt ...
                  + Sy_1 * sqrt(dt) * dW_x + Sy_2 * sqrt(dt) * dW_y;
    L_y(j) = - delta * y(j);
    f_y(j) = gamma * y(j);
    L_x(j) = beta * y(j) - alpha;
end

% Verify the claimed equivalence of the two branches' simulations rather than
% assuming it: re-run under the direction y -> x formulation and diff.
rng(42)
x2 = zeros(1, N+1); y2 = zeros(1, N+1);
x2(1) = 4; y2(1) = 4;
L_x2 = zeros(1, N+1); f_x2 = zeros(1, N+1); L_y2 = zeros(1, N+1); f_y2 = 0;
L_x2(1) = beta * x2(1);
f_x2(1) = - alpha * x2(1);
L_y2(1) = gamma - delta * x2(1);
for j = 2:N+1
    dW_x = randn;
    dW_y = randn;
    x2(j) = x2(j-1) + (L_x2(j-1) * y2(j-1) + f_x2(j-1)) * dt ...
                    + Sx_1 * sqrt(dt) * dW_x + Sx_2 * sqrt(dt) * dW_y;
    y2(j) = y2(j-1) + (L_y2(j-1) * y2(j-1) + f_y2) * dt ...
                    + Sy_1 * sqrt(dt) * dW_x + Sy_2 * sqrt(dt) * dW_y;
    L_x2(j) = beta * x2(j);
    f_x2(j) = - alpha * x2(j);
    L_y2(j) = gamma - delta * x2(j);
end
branch_gap = max([max(abs(x - x2)), max(abs(y - y2))]);
fprintf('BRANCH_EQUIVALENCE_GAP %.6e\n', branch_gap);

% ---- Generic scalar CGNS filter / smoother / metric -------------------------
% Written once with the roles passed in, exactly as the reference writes them
% twice with the names swapped. obs is the observed path; Lo couples the latent
% into the observed drift; fo is the rest of the observed drift; Ll is the
% latent self-drift; fl is the rest of the latent drift.
function [fm, fc, sm, sc, aci] = cgns(obs, Lo, fo, Ll, fl, Soo, Sll, mu_init, R_init, dt, N)
    Soo_inv = 1 / Soo;
    Sol = 0;  Slo = 0;   % the reference sets the noise cross-terms to zero here
    fm = zeros(1, N+1); fc = zeros(1, N+1);
    fm(1) = mu_init; fc(1) = R_init;
    mu0 = mu_init; R0 = R_init;
    for j = 2:N+1
        dobs = obs(j) - obs(j-1);
        aux = Slo + R0 * Lo(j-1);
        mu = mu0 + (Ll(j-1) * mu0 + fl(j-1)) * dt ...
                 + aux * Soo_inv * (dobs - (Lo(j-1) * mu0 + fo(j-1)) * dt);
        R = R0 + (Ll(j-1) * R0 + R0 * Ll(j-1) + Sll - aux * Soo_inv * aux) * dt;
        fm(j) = mu; fc(j) = R; mu0 = mu; R0 = R;
    end
    sm = zeros(1, N+1); sc = zeros(1, N+1);
    sm(N+1) = fm(N+1); sc(N+1) = fc(N+1);
    muT = sm(N+1); RT = sc(N+1);
    for j = N:-1:1
        dobs = obs(j+1) - obs(j);
        A_j = Ll(j) - Slo * Soo_inv * Lo(j);
        B_j = Sll - Slo * Soo_inv * Sol;
        mu = muT - (Ll(j) * muT + fl(j) - B_j / fc(j) * (fm(j) - muT)) * dt ...
                 + Slo * Soo_inv * (-dobs + (Lo(j) * muT + fo(j)) * dt);
        R = RT - ((A_j + B_j / fc(j)) * RT + RT * (A_j + B_j / fc(j)) - B_j) * dt;
        sm(j) = mu; sc(j) = R; muT = mu; RT = R;
    end
    ratio = sc ./ fc;
    aci = 0.5 * (sm - fm).^2 ./ fc + 0.5 * (-log(ratio) + ratio - 1);
end

idx = 1:40:(N+1);
zerovec = zeros(1, N+1);

% ---- Direction x -> y : prey observed, predator latent ----------------------
% Latent self-drift beta*y - alpha is time-varying: this is the path the dyad
% fixture cannot reach.
[fm1, fc1, sm1, sc1, aci1] = cgns( ...
    y, -delta * y, gamma * y, beta * y - alpha, zerovec, ...
    sigma_y^2, sigma_x^2, x(1), 0.1, dt, N);

% ---- Direction y -> x : predator observed, prey latent ----------------------
[fm2, fc2, sm2, sc2, aci2] = cgns( ...
    x, beta * x, -alpha * x, gamma - delta * x, zerovec, ...
    sigma_x^2, sigma_y^2, y(1), 0.1, dt, N);

% ---- Write fixtures ---------------------------------------------------------
T = array2table([(0:N)' * dt, x', y'], ...
                'VariableNames', {'t', 'predator', 'prey'});
writetable(T, 'predprey_signal.csv');

R1 = array2table([idx' , fm1(idx)', fc1(idx)', sm1(idx)', sc1(idx)', aci1(idx)'], ...
    'VariableNames', {'index','filter_mean','filter_cov','smoother_mean','smoother_cov','ACI_metric'});
writetable(R1, 'predprey_reference_predator_to_prey.csv');

R2 = array2table([idx', fm2(idx)', fc2(idx)', sm2(idx)', sc2(idx)', aci2(idx)'], ...
    'VariableNames', {'index','filter_mean','filter_cov','smoother_mean','smoother_cov','ACI_metric'});
writetable(R2, 'predprey_reference_prey_to_predator.csv');

fprintf('PREDPREY_ORACLE_DONE rows=%d Lyrange1=[%.4f %.4f] Lyrange2=[%.4f %.4f]\n', ...
        numel(idx), min(beta*y - alpha), max(beta*y - alpha), ...
        min(gamma - delta*x), max(gamma - delta*x));
