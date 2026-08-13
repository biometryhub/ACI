% ACI_ORACLE_CIR  Independent-oracle fixture for the causal influence range.
%
% Transcribes the online-smoother and causal-influence-range sections of the
% reference implementation (dyad_interaction_model.m, upstream commit
% 733c49fc5a36a5dad857608ef550a8f836466144) and writes the reference values
% aciR is graded against.
%
% The observed signal is NOT regenerated here. It is read from the dyad signal
% fixture already shipped in the package, whose bytes are hash-pinned in
% inst/extdata/oracle-manifest.yml, so this harness and the existing dyad
% oracle are driven by the identical path and no second signal fixture is
% introduced. Only the leading M steps are used: the reference online-smoother
% algorithm stores a staggered triangle that is quadratic in the record, which
% at the full 30001 steps would need several gigabytes to hold quantities that
% are immediately reduced to scalars.
%
% Run from the oracle/ directory:
%   matlab -nodisplay -nosplash -batch "run('aci_oracle_cir.m')"

clearvars; close all;

% simps.m is the reference implementation's quadrature helper and lives beside
% the upstream sources rather than here.
addpath('../matlab_reference');

% ---- Parameters of the reference dyad ---------------------------------------
d_x = 0.5; d_y = 0.5; gamma = 2; F_x = 0.5; F_y = 1;
sigma_x = 0.5; sigma_y = 1;
dt = 0.001;

S_xoS_x = sigma_x^2;
S_yoS_y = sigma_y^2;
S_yoS_x = 0;
S_xoS_y = 0;
S_xoS_x_inv = 1 / S_xoS_x;

% ---- Observed signal, read from the pinned fixture --------------------------
raw = readmatrix('../aciR/inst/extdata/dyad_signal_x.csv');
M = 2001;
x = raw(1:M, 2).';
N = M - 1;

L_x = gamma * x;
f_x = F_x - d_x * x;
L_y = -d_y;
f_y = F_y - gamma * x.^2;

% ---- Filter -----------------------------------------------------------------
filter_mean = zeros(1, N+1);
filter_cov = zeros(1, N+1);
filter_mean(1) = F_y / d_y;
filter_cov(1) = 0.1;
for j = 2:N+1
    dx = x(j) - x(j-1);
    aux = S_yoS_x + filter_cov(j-1) * L_x(j-1);
    filter_mean(j) = filter_mean(j-1) ...
        + (L_y * filter_mean(j-1) + f_y(j-1)) * dt ...
        + aux * S_xoS_x_inv * (dx - (L_x(j-1) * filter_mean(j-1) + f_x(j-1)) * dt);
    filter_cov(j) = filter_cov(j-1) ...
        + (2 * L_y * filter_cov(j-1) + S_yoS_y - aux * S_xoS_x_inv * aux) * dt;
end

% ---- Online-smoother auxiliary matrices -------------------------------------
E_j_matrices = zeros(1, N+1);
F_j_matrices = zeros(1, N+1);
for j = 1:N+1
    G_x_j = L_x(j) + S_xoS_y / filter_cov(j);
    G_y_j = L_y + S_yoS_y / filter_cov(j);
    C_jj = 1 - G_y_j * dt;
    H_j = filter_cov(j) \ (L_y * filter_cov(j) + filter_cov(j) * L_y + S_yoS_y);
    K_j = S_xoS_x_inv * G_x_j;
    E_j_matrices(j) = C_jj + S_yoS_x * K_j * dt;
    F_j_matrices(j) = - filter_cov(j) * ( ...
        K_j + (G_x_j * K_j * filter_cov(j) * K_j ...
               - filter_cov(j) \ H_j * filter_cov(j) * K_j + L_y * K_j) * dt ...
        - L_x(j) * (S_xoS_x_inv + K_j * filter_cov(j) * K_j * dt) );
end

% ---- Fixed-lag online smoother, full lag ------------------------------------
% Nested cell arrays exactly as the reference: row n holds the estimates at
% every time instant given observations up to n.
fixed_lag = N + 1;
online_fixed_mean = cell(N+1, 1);
online_fixed_cov = cell(N+1, 1);
update_matrices_fixed = cell(max(N-1, 1), 1);
for n = 1:(N-1)
    update_matrices_fixed{n} = zeros(1, n+1);
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end
for n = N:N+1
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end

online_fixed_mean{1}(1) = filter_mean(1);
online_fixed_cov{1}(1) = filter_cov(1);

online_fixed_mean{2}(2) = filter_mean(2);
online_fixed_cov{2}(2) = filter_cov(2);
aux_vec = filter_mean(1) ...
    - E_j_matrices(1) * ((1 + L_y * dt) * filter_mean(1) + f_y(1) * dt) ...
    + F_j_matrices(1) * (x(2) - x(1) - (L_x(1) * filter_mean(1) + f_x(1)) * dt);
online_fixed_mean{2}(1) = E_j_matrices(1) * filter_mean(2) + aux_vec;
aux_mat = filter_cov(1) ...
    - E_j_matrices(1) * (1 + L_y * dt) * filter_cov(1) ...
    - F_j_matrices(1) * L_x(1) * filter_cov(1) * dt;
online_fixed_cov{2}(1) = E_j_matrices(1) * filter_cov(2) * E_j_matrices(1) + aux_mat;

for n = 3:N+1
    online_fixed_mean{n}(n) = filter_mean(n);
    online_fixed_cov{n}(n) = filter_cov(n);

    aux_vec = filter_mean(n-1) ...
        - E_j_matrices(n-1) * ((1 + L_y * dt) * filter_mean(n-1) + f_y(n-1) * dt) ...
        + F_j_matrices(n-1) * (x(n) - x(n-1) - (L_x(n-1) * filter_mean(n-1) + f_x(n-1)) * dt);
    online_fixed_mean{n}(n-1) = E_j_matrices(n-1) * filter_mean(n) + aux_vec;
    aux_mat = filter_cov(n-1) ...
        - E_j_matrices(n-1) * (1 + L_y * dt) * filter_cov(n-1) ...
        - F_j_matrices(n-1) * L_x(n-1) * filter_cov(n-1) * dt;
    online_fixed_cov{n}(n-1) = E_j_matrices(n-1) * filter_cov(n) * E_j_matrices(n-1) + aux_mat;

    for j = (n-1):-1:1
        if (1 <= j) && (j <= n-1-fixed_lag)
            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j);
        elseif (n-fixed_lag <= j) && (j <= n-1)
            if j == n-1
                update_matrices_fixed{n-2}(n-1) = 1;
            elseif j == n-2
                update_matrices_fixed{n-2}(n-2) = E_j_matrices(n-2);
            else
                update_matrices_fixed{n-2}(j) = update_matrices_fixed{n-3}(j) * E_j_matrices(n-2);
            end
            online_mean_inov = online_fixed_mean{n}(n-1) - filter_mean(n-1);
            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j) + update_matrices_fixed{n-2}(j) * online_mean_inov;
            online_cov_inov = online_fixed_cov{n}(n-1) - filter_cov(n-1);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j) + update_matrices_fixed{n-2}(j) * online_cov_inov * update_matrices_fixed{n-2}(j);
        end
    end
end

% ---- Reference table: online smoother at sampled (j, n) pairs ---------------
j_samples = 51:100:1451;
rows = [];
for j = j_samples
    for n = [j, j+10, j+50, j+200, N+1]
        if n <= N+1
            rows(end+1, :) = [j, n, online_fixed_mean{n}(j), online_fixed_cov{n}(j)]; %#ok<SAGROW>
        end
    end
end
T1 = array2table(rows, 'VariableNames', {'j', 'n', 'online_mean', 'online_cov'});
writetable(T1, 'cir_online_reference.csv');

% ---- Reference table: divergence sequence and the ranges --------------------
eps_values = [1e-1, 1e-2, 1e-3, 1e-4];
out = zeros(numel(j_samples), 3 + numel(eps_values));
for idx = 1:numel(j_samples)
    j = j_samples(idx);
    RE_n = zeros(1, N+1-j+1);
    for obs = j:(N+1)
        cov_ratio = online_fixed_cov{end}(j) / online_fixed_cov{obs}(j);
        RE_n(obs-j+1) = 0.5 * (online_fixed_mean{end}(j) - online_fixed_mean{obs}(j))^2 / online_fixed_cov{obs}(j) ...
            + 0.5 * (cov_ratio - 1 - log(cov_ratio));
    end
    peak = max(RE_n);
    if peak > 1e-5
        objective = simps(RE_n) * dt / peak;
    else
        objective = 0;
    end
    subj = zeros(1, numel(eps_values));
    for e = 1:numel(eps_values)
        k = find(RE_n > eps_values(e), 1, 'last');
        if isempty(k)
            k = 0;
        end
        subj(e) = k * dt;
    end
    out(idx, :) = [j, peak, objective, subj];
end
names = [{'j', 'peak', 'objective'}, ...
         arrayfun(@(e) sprintf('subj_%g', e), eps_values, 'UniformOutput', false)];
T2 = array2table(out, 'VariableNames', names);
writetable(T2, 'cir_range_reference.csv');

fprintf('CIR_ORACLE_DONE M=%d rows1=%d rows2=%d\n', M, height(T1), height(T2));
