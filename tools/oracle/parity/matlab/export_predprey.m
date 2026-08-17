function export_predprey(workspace_path, direction, bundle_dir, report_dir)
%EXPORT_PREDPREY  Write a predator-prey bundle and its reference side.
%
%   The dyad path builds a dataset and drives the extracted kernels on it.
%   Predator-prey cannot be handled that way: `noisy_predator_prey_model.m`
%   must not be run top to bottom (see F6), so the graded reference values are
%   the ones the COMPOSED runner already produced and gate G1 already verified
%   the extracts reproduce exactly. This helper therefore reads the captured
%   workspace and writes two things: the model in the orientation aciR expects,
%   and the reference's own answers beside it.
%
%   The relabelling is the one place a mistake could hide, so it is stated
%   rather than buried. The reference names its variables after the ROLE each
%   plays in the causal direction under study, and the two directions are
%   mirror images:
%
%     direction 1, x(t) -> y : y (prey)     is OBSERVED, x (predator) is latent
%     direction 2, y(t) -> x : x (predator) is OBSERVED, y (prey)     is latent
%
%   aciR always calls the observed process x and the latent one y, which is the
%   dyad's convention and therefore direction two's. Direction one is the
%   mirror, and its columns are exchanged here (once, in MATLAB, next to the
%   discriminator that proves which direction was captured).
%
%   Nothing numerical is recomputed. The Grammians written below are the values
%   the reference's own filter block formed, not a re-derivation from the
%   feedback matrices, so a slot assigned to the wrong role shows up as a gross
%   disagreement in the filter rather than as a plausible number.
%
%   DIRECTION is 1 or 2. BUNDLE_DIR receives meta.dcf and arrays.csv;
%   REPORT_DIR receives the matlab_predprey_dir<D>*.csv reference side.

w = load(workspace_path);
name = sprintf('predprey_dir%d', direction);

% ---- the discriminator, before anything is written --------------------------
%
% `f_x` is a SCALAR in direction one (reference line 97) and an ARRAY in
% direction two (line 126); `f_y` is the mirror. This is the same condition the
% capture profile asserts, repeated here because this function performs a
% direction-dependent relabelling and must fail loudly rather than silently
% mirror the wrong workspace.
if direction == 1
    if numel(w.f_x) ~= 1 || numel(w.f_y) ~= w.N + 1
        error('parity:wrongDirection', ...
              ['workspace does not look like direction one: numel(f_x) = ' ...
               '%d, numel(f_y) = %d.'], numel(w.f_x), numel(w.f_y));
    end
    observed  = w.y;
    L_x_aci   = w.L_y;                   % latent -> observed coefficient
    f_x_aci   = w.f_y;                   % forcing in the observed process
    L_y_aci   = w.L_x;                   % latent self-drift
    f_y_aci   = w.f_x;                   % forcing in the latent process
    S_xoS_x   = w.S_yoS_y;               % observation-noise Grammian
    S_yoS_y   = w.S_xoS_x;               % latent-noise Grammian
    S_yoS_x   = w.S_yoS_x;
    S_xoS_y   = w.S_xoS_y;
elseif direction == 2
    if numel(w.f_y) ~= 1 || numel(w.f_x) ~= w.N + 1
        error('parity:wrongDirection', ...
              ['workspace does not look like direction two: numel(f_x) = ' ...
               '%d, numel(f_y) = %d.'], numel(w.f_x), numel(w.f_y));
    end
    observed  = w.x;
    L_x_aci   = w.L_x;
    f_x_aci   = w.f_x;
    L_y_aci   = w.L_y;
    f_y_aci   = w.f_y;
    S_xoS_x   = w.S_xoS_x;
    S_yoS_y   = w.S_yoS_y;
    S_yoS_x   = w.S_yoS_x;
    S_xoS_y   = w.S_xoS_y;
else
    error('parity:badDirection', 'direction must be 1 or 2, not %d.', ...
          direction);
end

N = w.N;
dt = w.dt;
% The scalar forcing of the latent process is broadcast, which is what the
% reference's own loops do with it implicitly.
if numel(f_y_aci) == 1
    f_y_aci = repmat(f_y_aci, 1, N + 1);
end
if numel(f_x_aci) == 1
    f_x_aci = repmat(f_x_aci, 1, N + 1);
end

% ---- the bundle -------------------------------------------------------------
%
% Grammians rather than the feedback matrices Sx_1..Sy_2 that the dyad bundle
% carries. Under the mirror those four names change role as well as value, and
% re-deriving a Grammian from them on the R side would be a second chance to
% get the mirror wrong; the reference computed these itself.
if ~isfolder(bundle_dir); mkdir(bundle_dir); end
writetable( ...
    table((0:N)' * dt, observed(:), L_x_aci(:), f_x_aci(:), L_y_aci(:), ...
          f_y_aci(:), ...
          'VariableNames', {'t', 'x', 'L_x', 'f_x', 'L_y', 'f_y'}), ...
    fullfile(bundle_dir, 'arrays.csv'));

plot_end_idx = round(w.time_end_plot / dt) + 1;

fid = fopen(fullfile(bundle_dir, 'meta.dcf'), 'w');
fprintf(fid, 'Name: %s\n', name);
fprintf(fid, 'Kind: scalar\n');
fprintf(fid, 'Direction: %d\n', direction);
if direction == 1
    fprintf(fid, ['Description: Predator-prey, causal direction x(t) -> y. ' ...
                  'The prey y is observed and the predator x is latent, so ' ...
                  'the reference columns are exchanged into aciR''s ' ...
                  'observed-is-x orientation.\n']);
else
    fprintf(fid, ['Description: Predator-prey, causal direction y(t) -> x. ' ...
                  'The predator x is observed and the prey y is latent, ' ...
                  'which is already aciR''s orientation.\n']);
end
fprintf(fid, ['Provenance: noisy_predator_prey_model.m, composed runner ' ...
              'for profile %s, captured workspace\n'], name);
fprintf(fid, 'N: %d\n', N);
fprintf(fid, 'dt: %.17g\n', dt);
fprintf(fid, 'mu0: %.17g\n', w.filter_mean(1));
fprintf(fid, 'R0: %.17g\n', w.filter_cov(1));
fprintf(fid, 'S_xoS_x: %.17g\n', S_xoS_x);
fprintf(fid, 'S_yoS_y: %.17g\n', S_yoS_y);
fprintf(fid, 'S_yoS_x: %.17g\n', S_yoS_x);
fprintf(fid, 'S_xoS_y: %.17g\n', S_xoS_y);
fprintf(fid, 'CIRStart: %.17g\n', w.time_start_plot);
fprintf(fid, 'CIREnd: %.17g\n', w.time_end_plot);
fprintf(fid, 'FirstIdx: %d\n', w.first_idx);
fprintf(fid, 'LastIdx: %d\n', w.last_idx);
fprintf(fid, 'PlotLen: %d\n', w.plot_len);
% The end of the region the reference PLOTS, as distinct from the end of the
% region it computes. Exported rather than re-derived on the R side so the
% branch that produced last_idx is read off the run instead of retyped.
fprintf(fid, 'PlotEndIdx: %d\n', plot_end_idx);
fprintf(fid, 'Lookahead: %.17g\n', w.lookahead_tolerance);
fprintf(fid, 'EpsilonResolution: %d\n', w.epsilon_resolution);
fprintf(fid, 'FixedLag: %d\n', w.fixed_lag);
fprintf(fid, 'DefnLen: %d\n', numel(w.defn_objective_CIR));
fclose(fid);

% ---- the reference side -----------------------------------------------------
if ~isfolder(report_dir); mkdir(report_dir); end
stem = fullfile(report_dir, sprintf('matlab_%s', name));

writetable( ...
    table((0:N)' * dt, observed(:), w.filter_mean(:), w.filter_cov(:), ...
          w.smoother_mean(:), w.smoother_cov(:), w.ACI_metric(:), ...
          w.signal_smoother_filter(:), w.dispersion_smoother_filter(:), ...
          w.E_j_matrices(:), w.F_j_matrices(:), ...
          'VariableNames', {'t', 'x', 'filter_mean', 'filter_cov', ...
                            'smoother_mean', 'smoother_cov', 'ACI_metric', ...
                            'signal_part', 'dispersion_part', 'E_j', ...
                            'F_j'}), ...
    sprintf('%s.csv', stem));

% The reference runs the online smoother at one lag only (it sets
% fixed_lag = N+1 unconditionally), so that is the only lag there is anything
% to compare against, and the file is named with it rather than with a bare
% "full" so the comparison cannot drift onto a different quantity.
final_mean = w.online_fixed_mean{end}(:);
final_cov = w.online_fixed_cov{end}(:);
writetable( ...
    table((0:N)' * dt, final_mean, final_cov, ...
          'VariableNames', {'t', 'online_mean', 'online_cov'}), ...
    sprintf('%s_online_lag%d.csv', stem, w.fixed_lag));

writetable( ...
    table((w.first_idx:w.last_idx)', ((w.first_idx:w.last_idx)' - 1) * dt, ...
          w.approx_objective_CIR(:), w.max_RE_metric(:), ...
          'VariableNames', {'index', 't', 'objective', 'peak'}), ...
    sprintf('%s_cir.csv', stem));
writematrix(w.subjective_CIR, sprintf('%s_cir_subjective.csv', stem));
writematrix(10 .^ w.eps_ord_values(:), sprintf('%s_cir_epsilon.csv', stem));
writematrix(w.defn_objective_CIR(:), sprintf('%s_cir_defn.csv', stem));

fprintf(['exported %s: N = %d, window %d:%d (plot_len %d), plotted to %d, ' ...
         '%d thresholds, defn over %d times\n'], ...
        name, N, w.first_idx, w.last_idx, w.plot_len, plot_end_idx, ...
        w.epsilon_resolution, numel(w.defn_objective_CIR));

end
