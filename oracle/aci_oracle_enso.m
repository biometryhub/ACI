% ACI_ORACLE_ENSO  Independent-oracle fixture for the stochastic ENSO model.
%
% Transcribes the coefficient construction and the vector CGNS filter/smoother
% of ENSO_model_cond_ACI_u_h_W_tau_unobs.m (upstream commit
% 733c49fc5a36a5dad857608ef550a8f836466144), for the configuration in which
% (T_C, T_E, I) are observed and (u, h_W, tau) are not.
%
% The observed paths are generated here by a plain Euler-Maruyama integration
% rather than the reference's mixed scheme. That is deliberate and it is the
% one place this harness departs from upstream: the reference's wind-burst
% update carries a Milstein correction using the derivative of that variable's
% diffusion with respect to a DIFFERENT variable, which is not the standard
% scheme and which the accompanying paper never mentions. See
% design/2026-08-13_milstein_anomaly.md. The observed path here is therefore a
% path, not THE reference path, which is all an oracle for the FILTER needs,
% since the filter is graded on whatever signal it is driven with.
%
% What this fixture grades is the coefficient construction and the vector
% recursions on a system whose coefficients are matrix-valued, seasonally
% modulated, and whose BOTH noise covariances vary in time.
%
% Run from the oracle/ directory:
%   matlab -nodisplay -nosplash -batch "run('aci_oracle_enso.m')"

clearvars; close all;
rng(1618)

N = 4000;
dt = 0.005;

% ---- Parameters, exactly as the reference derives them ----------------------
factor = 0.65;
b_0 = 2.5; mu = 0.5;
alpha_2 = 0.125 * factor;
alpha_1 = alpha_2/2 * factor;
delta_u = alpha_1 * b_0 * mu;
delta_h = alpha_2 * b_0 * mu;
r = 0.25 * factor;
gamma_C = 0.75 * factor;
gamma_E = 0.75 * factor;
r_C = gamma_C*b_0*mu/2;
r_E = 3*gamma_E*b_0*mu/2;
zeta_C = gamma_C*b_0*mu/2;
zeta_E = gamma_E*b_0*mu/2;
C_u = 0.03 * factor;
d_tau = 2;
lambda = 2/60;
m = 2;
sigma_u = 0.04 * sqrt(factor);
sigma_h = 0.02 * sqrt(factor);
sigma_C = 0.04 * sqrt(factor);
sigma_E = sqrt(5) * 1e-2 * sqrt(factor);

% ---- Simulate the six-dimensional system ------------------------------------
u = zeros(1,N+1); h_W = zeros(1,N+1); T_C = zeros(1,N+1);
T_E = zeros(1,N+1); tau = zeros(1,N+1); I = zeros(1,N+1);
u(1) = 6.9136e-04; h_W(1) = -0.0028; T_C(1) = 0.0039;
T_E(1) = 0.0051; tau(1) = -0.0256; I(1) = 1.5841;

for j = 2:N+1
    t = (j-2)*dt;
    dW = randn(6,1) * sqrt(dt);

    sigma_I = sqrt(lambda * (4-I(j-1)) * I(j-1));
    I(j) = I(j-1) - lambda*(I(j-1) - m)*dt + sigma_I*dW(1);
    % Keep the Walker circulation strictly inside its domain: the noise
    % vanishes at both ends and the filter inverts a covariance built from it.
    I(j) = min(max(I(j), 0.05), 3.95);

    sigma = I(j-1)/5 * factor;
    c_1 = (25*(T_C(j-1) + 0.75/7.5).^2 + 0.9) * (1 + 0.3*sin(t*2*pi/6 - pi/6)) * factor;
    c_2 = 1.4*factor*(1 + 0.3*sin(t*2*pi/6 + 2*pi/6) + 0.25*sin(2*t*2*pi/6 + 2*pi/6));
    burst = (1 + (1 - I(j-1)/5))*0.15*sqrt(factor);
    beta_u = -0.2*burst; beta_h = -0.4*burst;
    beta_C = 0.8*burst;  beta_E = 1*burst;

    u(j)   = u(j-1)   + (-r*u(j-1)   - delta_u*(T_C(j-1)+T_E(j-1))/2)*dt + beta_u*tau(j-1)*dt + sigma_u*dW(2);
    h_W(j) = h_W(j-1) + (-r*h_W(j-1) - delta_h*(T_C(j-1)+T_E(j-1))/2)*dt + beta_h*tau(j-1)*dt + sigma_h*dW(3);
    T_C(j) = T_C(j-1) + ((r_C - c_1)*T_C(j-1) + zeta_C*T_E(j-1) + gamma_C*h_W(j-1) + sigma*u(j-1) + C_u)*dt + beta_C*tau(j-1)*dt + sigma_C*dW(4);
    T_E(j) = T_E(j-1) + ((r_E - c_2)*T_E(j-1) - zeta_E*T_C(j-1) + gamma_E*h_W(j-1))*dt + beta_E*tau(j-1)*dt + sigma_E*dW(5);

    sigma_tau = 0.9*(tanh(7.5*T_C(j-1)) + 1)*(1 + 0.3*cos(t*2*pi/6 + 2*pi/6));
    tau(j) = tau(j-1) - d_tau*tau(j-1)*dt + sigma_tau*dW(6);
end

% ---- Coefficients, as the reference builds them -----------------------------
L_x = zeros(3,3,N+1); L_y = zeros(3,3,N+1);
f_x = zeros(3,N+1);   f_y = zeros(3,N+1);
S_xoS_x = zeros(3,3,N+1); S_yoS_y = zeros(3,3,N+1);

for j = 1:N+1
    t = (j-1)*dt;
    burst = (1 + (1 - I(j)/5))*0.15*sqrt(factor);
    c_1 = (25*(T_C(j) + 0.75/7.5).^2 + 0.9) * (1 + 0.3*sin(t*2*pi/6 - pi/6)) * factor;
    c_2 = 1.4*factor*(1 + 0.3*sin(t*2*pi/6 + 2*pi/6) + 0.25*sin(2*t*2*pi/6 + 2*pi/6));

    L_x(:,:,j) = [ I(j)/5*factor, gamma_C, 0.8*burst;
                   0,             gamma_E, 1*burst;
                   0,             0,       0 ];
    f_x(:,j) = [ (r_C - c_1)*T_C(j) + zeta_C*T_E(j) + C_u;
                 (r_E - c_2)*T_E(j) - zeta_E*T_C(j);
                 -lambda*(I(j) - m) ];
    L_y(:,:,j) = [ -r, 0,  -0.2*burst;
                    0, -r, -0.4*burst;
                    0,  0, -d_tau ];
    f_y(:,j) = [ -delta_u*(T_C(j)+T_E(j))/2;
                 -delta_h*(T_C(j)+T_E(j))/2;
                 0 ];
    S_xoS_x(:,:,j) = diag([sigma_C^2, sigma_E^2, lambda*(4-I(j))*I(j)]);
    sigma_tau = 0.9*(tanh(7.5*T_C(j)) + 1)*(1 + 0.3*cos(t*2*pi/6 + 2*pi/6));
    S_yoS_y(:,:,j) = diag([sigma_u^2, sigma_h^2, sigma_tau^2]);
end

% ---- Filter and smoother, vector CGNS ---------------------------------------
x = [T_C; T_E; I];
fm = zeros(3,N+1); fc = zeros(3,3,N+1);
fm(:,1) = [u(1); h_W(1); tau(1)];
fc(:,:,1) = 0.01 * eye(3);
mu_f = fm(:,1); R = fc(:,:,1);

for j = 2:N+1
    inv_S = inv(S_xoS_x(:,:,j-1));
    dx = x(:,j) - x(:,j-1);
    aux = R * L_x(:,:,j-1)';
    mu_f = mu_f + (L_y(:,:,j-1)*mu_f + f_y(:,j-1))*dt ...
           + aux*inv_S*(dx - (L_x(:,:,j-1)*mu_f + f_x(:,j-1))*dt);
    R = R + (L_y(:,:,j-1)*R + R*L_y(:,:,j-1)' + S_yoS_y(:,:,j-1) - aux*inv_S*aux')*dt;
    R = (R + R')/2;
    fm(:,j) = mu_f; fc(:,:,j) = R;
end

sm = zeros(3,N+1); sc = zeros(3,3,N+1);
sm(:,N+1) = fm(:,N+1); sc(:,:,N+1) = fc(:,:,N+1);
muT = sm(:,N+1); RT = sc(:,:,N+1);
for j = N:-1:1
    A_j = L_y(:,:,j);
    B_j = S_yoS_y(:,:,j);
    BRf = B_j / fc(:,:,j);
    muT = muT - (L_y(:,:,j)*muT + f_y(:,j) - BRf*(fm(:,j) - muT))*dt;
    M = A_j + BRf;
    RT = RT - (M*RT + RT*M' - B_j)*dt;
    RT = (RT + RT')/2;
    sm(:,j) = muT; sc(:,:,j) = RT;
end

aci = zeros(1,N+1);
for j = 1:N+1
    Rf = fc(:,:,j); Rs = sc(:,:,j);
    d = sm(:,j) - fm(:,j);
    aci(j) = 0.5*(d'/Rf*d + trace(Rf\Rs) - 3 - log(det(Rs)/det(Rf)));
end

% ---- Write fixtures ---------------------------------------------------------
idx = 1:20:(N+1);
signal = array2table([(0:N)'*dt, T_C', T_E', I'], ...
                     'VariableNames', {'t','T_C','T_E','I'});
writetable(signal, 'enso_signal.csv');

ref = array2table([idx', fm(1,idx)', fm(2,idx)', fm(3,idx)', ...
    squeeze(fc(1,1,idx)), squeeze(fc(3,3,idx)), squeeze(fc(1,3,idx)), ...
    sm(1,idx)', sm(3,idx)', squeeze(sc(1,1,idx)), squeeze(sc(3,3,idx)), aci(idx)'], ...
    'VariableNames', {'index','fm1','fm2','fm3','fc11','fc33','fc13', ...
                      'sm1','sm3','sc11','sc33','ACI_metric'});
writetable(ref, 'enso_reference.csv');

fprintf('ENSO_ORACLE_DONE rows=%d Irange=[%.3f %.3f] Svar=[%.3e %.3e] metric=[%.4f %.4f]\n', ...
        numel(idx), min(I), max(I), min(squeeze(S_xoS_x(3,3,:))), ...
        max(squeeze(S_xoS_x(3,3,:))), min(aci), max(aci));
