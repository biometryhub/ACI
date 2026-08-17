function export_vector_bundle(workspace_path, bundle_dir, M)
%EXPORT_VECTOR_BUNDLE  Write a vector dataset bundle from a captured workspace.
%
%   The vector counterpart of the scalar contract, and it faces one question the
%   scalar case does not: how to put a 3-D array of matrices into CSV. Three
%   options were considered.
%
%     (a) One row per time step, the slice flattened column-major (chosen).
%         MATLAB's A(:,:,j)(:) and R's as.vector(A[,,j]) are both column-major,
%         so neither side transposes and a transposition bug has nowhere to
%         hide. n rows by r*c columns.
%     (b) Long format, one row per (slice, i, j, value). Self-describing, but
%         three times the bytes and both sides must reassemble, which is more
%         code on the axis where the mistakes actually happen.
%     (c) One file per matrix entry, A11.csv .. A33.csv. Trivially readable but
%         nine files per component, and the naming becomes the contract.
%
%   The record is truncated to M steps. That is sound for grading because BOTH
%   sides then run on the same M-step record: the comparison is aciR against
%   MATLAB, not against the published figure.
%
%   Noise is carried as the feedback matrices S_x and S_y, not as the
%   Grammians, matching the scalar bundle, and deliberately, because forming
%   the Grammians and their inverse is itself a place the two implementations
%   diverge (the reference pseudo-inverts, aciR factorises).

w = load(workspace_path);
if ~isfolder(bundle_dir); mkdir(bundle_dir); end
flat = @(A) reshape(A(:, :, 1:M), [], M).';        % M rows, r*c cols, col-major

writematrix([w.T_C(1:M)', w.T_E(1:M)', w.I(1:M)'], fullfile(bundle_dir, 'x.csv'));
writematrix(flat(w.L_x), fullfile(bundle_dir, 'L_x.csv'));
writematrix(flat(w.L_y), fullfile(bundle_dir, 'L_y.csv'));
writematrix(w.f_x(:, 1:M).', fullfile(bundle_dir, 'f_x.csv'));
writematrix(w.f_y(:, 1:M).', fullfile(bundle_dir, 'f_y.csv'));
writematrix(flat(w.S_x), fullfile(bundle_dir, 'S_x.csv'));
writematrix(flat(w.S_y), fullfile(bundle_dir, 'S_y.csv'));
writematrix([w.u(1); w.h_W(1); w.tau(1)].', fullfile(bundle_dir, 'mu0.csv'));

fid = fopen(fullfile(bundle_dir, 'meta.dcf'), 'w');
fprintf(fid, 'Name: enso_uhwtau_head\n');
fprintf(fid, 'Kind: vector\n');
fprintf(fid, 'Description: First %d steps of the reference ENSO run, (T_C,T_E,I) observed and (u,h_W,tau) unobserved.\n', M);
fprintf(fid, 'Provenance: ENSO_model_cond_ACI_u_h_W_tau_unobs.m lines 880-1163, unmodified\n');
fprintf(fid, 'N: %d\n', M - 1);
fprintf(fid, 'dt: %.17g\n', w.dt);
fprintf(fid, 'k: 3\n');
fprintf(fid, 'l: 3\n');
fprintf(fid, 'R0: 0.1\n');
% Whether S_x and S_y act on the SAME Wiener increments or on disjoint ones is
% a property of the model, not of the file format, and it decides whether the
% noise cross-Grammian is S_y*S_x' or exactly zero. It is declared rather than
% inferred: reading these two 3x3 blocks as sharing three increments gives a
% cross-Grammian of 0.38 here and destabilises the filter, where this model's
% own script (line 1231) states the cross terms are absent.
fprintf(fid, 'NoiseCoupling: disjoint\n');
fclose(fid);
fprintf('bundle written: %s (%d steps)\n', bundle_dir, M);

% ---- reference side, on the truncated record --------------------------------
N = M - 1; dt = w.dt;
T_C = w.T_C(1:M); T_E = w.T_E(1:M); I = w.I(1:M);
L_x = w.L_x(:, :, 1:M); L_y = w.L_y(:, :, 1:M);
f_x = w.f_x(:, 1:M); f_y = w.f_y(:, 1:M);
S_x = w.S_x(:, :, 1:M); S_y = w.S_y(:, :, 1:M);

[filter_mean, filter_cov, S_xoS_x_inv, S_yoS_y] = ref_enso_filter( ...
    T_C, T_E, I, w.u(1), w.h_W(1), w.tau(1), L_x, L_y, f_x, f_y, S_x, S_y, N, dt);
[smoother_mean, smoother_cov, E_j_matrices, F_j_matrices] = ref_enso_smoother( ...
    T_C, T_E, I, filter_mean, filter_cov, L_x, L_y, f_x, f_y, S_yoS_y, ...
    S_xoS_x_inv, N, dt);
[sig_p, dsp_p, ACI_metric] = ref_enso_metric(filter_mean, filter_cov, ...
                                             smoother_mean, smoother_cov);

out = fileparts(bundle_dir);
writematrix(filter_mean.',    fullfile(out, 'ref_enso_filter_mean.csv'));
writematrix(reshape(filter_cov, [], M).',   fullfile(out, 'ref_enso_filter_cov.csv'));
writematrix(smoother_mean.',  fullfile(out, 'ref_enso_smoother_mean.csv'));
writematrix(reshape(smoother_cov, [], M).', fullfile(out, 'ref_enso_smoother_cov.csv'));
writematrix(ACI_metric.',     fullfile(out, 'ref_enso_metric.csv'));

% ---- F8: the two routes to the observation-noise inverse --------------------
%
% The reference pseudo-inverts; aciR factorises. Both are computed here on the
% SAME Grammians so the difference is reported rather than assumed negligible.
worst_diff = 0; worst_cond = 0; n_singular = 0;
for j = 1:M
    G = S_x(:, :, j) * S_x(:, :, j)';
    P = pinv(G);
    worst_cond = max(worst_cond, cond(G));
    [R, flag] = chol(G);
    if flag == 0
        C = R \ (R' \ eye(3));
        worst_diff = max(worst_diff, max(abs(P(:) - C(:))));
    else
        n_singular = n_singular + 1;
    end
end
fprintf('F8 pinv vs chol: worst |diff| = %.6g over %d slices\n', worst_diff, M);
fprintf('F8 worst condition number = %.6g ; slices chol refused = %d\n', ...
        worst_cond, n_singular);

end
