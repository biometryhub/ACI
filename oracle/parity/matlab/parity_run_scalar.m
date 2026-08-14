function parity_run_scalar(dataset_dir, out_path, do_cir)
%PARITY_RUN_SCALAR  Run the extracted reference kernels on a dataset bundle.
%
%   The MATLAB half of the side-by-side. Reads the same bundle aciR reads,
%   drives the byte-verified extracts, and writes one CSV of per-step
%   quantities for R to compare against its own.
%
%   Two things about the interface are load-bearing and are asserted here
%   rather than assumed:
%
%   * The reference's filter reads the true latent path only as y(1), its
%     initial mean -- asserted by the extractor's ConsumedOnly rule. That is
%     what makes it legitimate to pass the scalar mu0 in its place on a dataset
%     that has no latent path at all.
%
%   * The reference hardcodes the initial covariance as 0.1 inside the
%     extracted range, so it is not an argument. A dataset declaring any other
%     R0 cannot be compared against this side, and is rejected rather than
%     silently graded at the wrong value.

if nargin < 3
    do_cir = true;
end

meta = read_meta(fullfile(dataset_dir, 'meta.dcf'));
arrays = readtable(fullfile(dataset_dir, 'arrays.csv'));

if abs(meta.R0 - 0.1) > 0
    error('parity:R0Mismatch', ...
          ['dataset declares R0 = %.17g, but the reference hardcodes 0.1 ' ...
           'inside ref_filter_scalar (dyad_interaction_model.m:150). ' ...
           'Comparison refused.'], meta.R0);
end

N = meta.N; dt = meta.dt;
x = arrays.x(:)'; L_x = arrays.L_x(:)'; f_x = arrays.f_x(:)';
f_y = arrays.f_y(:)'; L_y = meta.L_y;
y = meta.mu0;                       % consumed only as y(1); see above

[filter_mean, filter_cov, S_xoS_x, S_xoS_x_inv, S_yoS_y, S_yoS_x, S_xoS_y] = ...
    ref_filter_scalar(x, y, L_x, f_x, f_y, L_y, meta.Sx_1, meta.Sx_2, ...
                      meta.Sy_1, meta.Sy_2, N, dt);

[smoother_mean, smoother_cov, E_j_matrices, F_j_matrices] = ...
    ref_smoother_scalar(x, filter_mean, filter_cov, L_x, f_x, f_y, L_y, ...
                        S_xoS_x_inv, S_yoS_y, S_yoS_x, S_xoS_y, N, dt);

[ACI_metric, signal_part, cov_ratio, dispersion_part] = ...
    ref_aci_metric_scalar(filter_mean, filter_cov, smoother_mean, smoother_cov);

t = (0:N)' * dt;
core = table(t, x(:), filter_mean(:), filter_cov(:), smoother_mean(:), ...
             smoother_cov(:), ACI_metric(:), signal_part(:), ...
             dispersion_part(:), E_j_matrices(:), F_j_matrices(:), ...
    'VariableNames', {'t','x','filter_mean','filter_cov','smoother_mean', ...
                      'smoother_cov','ACI_metric','signal_part', ...
                      'dispersion_part','E_j','F_j'});
writetable(core, out_path);
fprintf('parity core: %d steps -> %s\n', N + 1, out_path);

% ---- online smoother and CIR ------------------------------------------------
%
% Quadratic in N on both time and memory, so they run only when asked. The
% online smoother is driven at three lags: zero, where it must collapse onto
% the filter, a short lag, and the full path, where it must reach the backward
% smoother.
if do_cir
    [~, cir_stem, ~] = fileparts(out_path);
    cir_dir = fileparts(out_path);
    lags = [0, 25, N + 1];
    for k = 1:numel(lags)
        fixed_lag = lags(k);
        [online_fixed_mean, online_fixed_cov] = ...
            ref_online_smoother_scalar(x, L_x, f_x, f_y, L_y, filter_mean, ...
                                       filter_cov, E_j_matrices, ...
                                       F_j_matrices, fixed_lag, N, dt);
        final_mean = online_fixed_mean{end}(:);
        final_cov = online_fixed_cov{end}(:);
        writetable( ...
            table((0:N)'*dt, final_mean, final_cov, ...
                  'VariableNames', {'t','online_mean','online_cov'}), ...
            fullfile(cir_dir, sprintf('%s_online_lag%d.csv', cir_stem, ...
                                      fixed_lag)));
        fprintf('parity online: lag %d -> %s_online_lag%d.csv\n', ...
                fixed_lag, cir_stem, fixed_lag);

        % The causal influence range is built on the FULL-lag online smoother,
        % which is the only lag the reference ever uses (it sets fixed_lag =
        % N+1 unconditionally). Running it at the other lags would compare a
        % quantity the reference does not define.
        if fixed_lag == N + 1
            T = N * dt;
            time_start_plot = meta.CIRStart;
            time_end_plot = meta.CIREnd;
            [subjective_CIR, approx_objective_CIR, RE_metric, max_RE_metric, ...
             eps_ord_values, epsilon_resolution, lookahead_tolerance, ...
             first_idx, last_idx, plot_len] = ...
                ref_cir_scalar(dt, T, time_start_plot, time_end_plot, ...
                               online_fixed_mean, online_fixed_cov);
            writetable( ...
                table((first_idx:last_idx)', ...
                      ((first_idx:last_idx)'-1)*dt, ...
                      approx_objective_CIR(:), max_RE_metric(:), ...
                      'VariableNames', {'index','t','objective','peak'}), ...
                fullfile(cir_dir, sprintf('%s_cir.csv', cir_stem)));
            writematrix(subjective_CIR, ...
                fullfile(cir_dir, sprintf('%s_cir_subjective.csv', cir_stem)));
            writematrix(10.^eps_ord_values(:), ...
                fullfile(cir_dir, sprintf('%s_cir_epsilon.csv', cir_stem)));
            fprintf(['parity cir: %d reporting times, %d thresholds, ' ...
                     'RE_metric %dx%d\n'], plot_len, epsilon_resolution, ...
                    size(RE_metric, 1), size(RE_metric, 2));
            clear RE_metric subjective_CIR
        end

        % The staggered triangles are quadratic in N and MATLAB does not
        % reclaim them promptly across loop iterations. Holding three lags'
        % worth at once killed the process at N = 3000; they are released as
        % soon as the quantities derived from them have been written.
        clear online_fixed_mean online_fixed_cov
    end
end

end
