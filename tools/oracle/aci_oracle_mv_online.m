% ACI_ORACLE_MV_ONLINE  Vector ONLINE SMOOTHER with matrix cross-noise.
%
% THIS HARNESS HAS NO UPSTREAM COUNTERPART, AND THE GAP IS LARGER HERE THAN
% ANYWHERE ELSE IN THE PACKAGE.
%
% Every scalar model in the reference implementation sets the noise
% cross-covariance to zero, and every ENSO script says so in as many words,
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
fm = zeros(n_y, N+1); fc = zeros(n_y, n_y, N+1);
fm(:,1) = [0.8; 0.2];
fc(:,:,1) = 0.2*eye(n_y);
mu = fm(:,1); R = fc(:,:,1);
for j = 2:N+1
    dx = x(:,j) - x(:,j-1);
    aux = S_yoS_x + R*L_x(:,:,j-1)';
    mu = mu + (L_y(:,:,j-1)*mu + f_y(:,j-1))*dt ...
            + aux*S_xoS_x_inv*(dx - (L_x(:,:,j-1)*mu + f_x(:,j-1))*dt);
    R = R + (L_y(:,:,j-1)*R + R*L_y(:,:,j-1)' + S_yoS_y - aux*S_xoS_x_inv*aux')*dt;
    R = (R+R')/2;
    fm(:,j) = mu; fc(:,:,j) = R;
end

% ---- Auxiliary matrices, equations (3.5)-(3.7) in FULL generality ------------
% The reference implementation carries only the zero-cross-noise
% specialisation, equation (3.8). This is an independent transcription of the
% general form; it is NOT an authors'-reference grounding.
S_xoS_y = S_yoS_x';
E = zeros(n_y, n_y, N+1); F = zeros(n_y, n_x, N+1);
for j = 1:N+1
    Rf = fc(:,:,j); Rfi = inv(Rf);
    G_x = L_x(:,:,j) + S_xoS_y*Rfi;                       % k x l
    G_y = L_y(:,:,j) + S_yoS_y*Rfi;                       % l x l
    H   = Rfi*(L_y(:,:,j)*Rf + Rf*L_y(:,:,j)' + S_yoS_y); % l x l
    K   = S_xoS_x_inv*G_x;                                % k x l
    E(:,:,j) = eye(n_y) + (S_yoS_x*S_xoS_x_inv*G_x - G_y)*dt;
    F(:,:,j) = -Rf * ( K' ...
        + (G_x'*K*Rf*K' - Rfi*H'*Rf*K' + L_y(:,:,j)'*K')*dt ...
        - L_x(:,:,j)'*(S_xoS_x_inv + K*Rf*K'*dt) );
end

% ---- Online smoother, equations (3.10)-(3.16) -------------------------------
innov_m = zeros(n_y, N); innov_c = zeros(n_y, n_y, N);
for k = 1:N
    Rf = fc(:,:,k);
    b = fm(:,k) - E(:,:,k)*((eye(n_y) + L_y(:,:,k)*dt)*fm(:,k) + f_y(:,k)*dt) ...
        + F(:,:,k)*(x(:,k+1) - x(:,k) - (L_x(:,:,k)*fm(:,k) + f_x(:,k))*dt);
    P = Rf - E(:,:,k)*(eye(n_y) + L_y(:,:,k)*dt)*Rf - F(:,:,k)*L_x(:,:,k)*Rf*dt;
    innov_m(:,k)   = E(:,:,k)*fm(:,k+1) + b - fm(:,k);
    innov_c(:,:,k) = E(:,:,k)*fc(:,:,k+1)*E(:,:,k)' + P - Rf;
end

j_samples = 51:200:1651;
rows = [];
for j = j_samples
    mu_j = fm(:,j); R_j = fc(:,:,j); D = eye(n_y);
    for n_obs = j:(N)
        if any(n_obs == [j, j+5, j+40, j+200, N])
            rows(end+1,:) = [j, n_obs, mu_j', R_j(1,1), R_j(1,2), R_j(2,2)]; %#ok<SAGROW>
        end
        mu_j = mu_j + D*innov_m(:,n_obs);
        R_j  = R_j  + D*innov_c(:,:,n_obs)*D';
        R_j  = (R_j + R_j')/2;
        D = D*E(:,:,n_obs);
    end
end
T = array2table(rows, 'VariableNames', {'j','n','om1','om2','oc11','oc12','oc22'});
writetable(T, 'mv_online_reference.csv');
fprintf('MV_ONLINE_DONE rows=%d crossnorm=%.4e\n', height(T), norm(S_yoS_x,'fro'));
