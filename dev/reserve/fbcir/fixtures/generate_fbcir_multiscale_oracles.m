% acir reserve file
% Origin: aci/tests/testthat/fixtures/oracles/generate_fbcir_multiscale_oracles.m
% Source package: aci 0.0.30, git tree 97f6b124
% Category: fbcir
% Intended release: 0.2.0 or 0.3.0 (FBCIR_code family; order TBD)
% Reason: MATLAB generator for the FBCIR multiscale fixtures; never run, referenced by no test.
% Verbatim copy from the aci 0.0.30 sources; not modified.

function generate_fbcir_multiscale_oracles(source_dir, output_dir, n_steps)
%GENERATE_FBCIR_MULTISCALE_ORACLES Execute the supplied FBCIR MATLAB source.
%
% This development-only harness creates compact numerical fixtures for all
% three multiscale routes in the authors' FBCIR repository:
%
%   (y1,y2) -> (x1,x2)
%   y1 -> x1 | (x2,y2)
%   y2 -> x2 | (x1,y1)
%
% It does not restate the filter, smoother, online-smoother, or CIR equations.
% Instead it verifies the exact SHA-256 of each supplied script, makes a
% temporary copy, and applies only these auditable development changes:
%
%   * shorten N (default 300 rather than 20000);
%   * start the reporting window at t = 0;
%   * hide figures and disable optional heatmap caches;
%   * retain online histories long enough to export sampled values.
%
% The numerical algorithms, seed, dt, model parameters, conditional precision
% masks, covariance regularisation, epsilon grid, quadrature, and indexing are
% otherwise executed by the supplied scripts themselves. Consequently the
% generated files are AUTHORS-SOURCE EXECUTION fixtures. They are not
% author-provided outputs and must never be described that way.
%
% Example (from this file's directory):
%
%   generate_fbcir_multiscale_oracles( ...
%       '../../../../../FBCIR_code-main', pwd, 300)
%
% Written for MATLAB batch execution. The supplied repository reports R2024b
% as its reference environment; the generation record captures the release
% actually used so a different release cannot be silently conflated with it.

if nargin < 3 || isempty(n_steps)
    n_steps = 300;
end
if nargin < 2 || isempty(output_dir)
    output_dir = pwd;
end
if nargin < 1 || isempty(source_dir)
    error('ACI:FBCIROracleSource', 'source_dir is required.');
end
validateattributes(n_steps, {'numeric'}, ...
    {'scalar', 'integer', '>=', 20, '<=', 1000});
source_dir = char(source_dir);
output_dir = char(output_dir);
if ~isfolder(source_dir)
    error('ACI:FBCIROracleSource', 'FBCIR source directory does not exist: %s', source_dir);
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end
[ok_source, source_attrs] = fileattrib(source_dir);
[ok_output, output_attrs] = fileattrib(output_dir);
if ~ok_source || ~ok_output
    error('ACI:FBCIROraclePath', 'Could not resolve source/output directories.');
end
source_dir = source_attrs.Name;
output_dir = output_attrs.Name;

routes = struct( ...
    'id', {'joint', 'y1', 'y2'}, ...
    'file', { ...
        'multiscale_joint_y_1_y_2_cause.m', ...
        'multiscale_marginal_y_1_cause.m', ...
        'multiscale_marginal_y_2_cause.m'}, ...
    'sha256', { ...
        '06685f0d46868bac65a47ff161a7b8c9933a2b0f1fb1a0f53132bd97e6b51827', ...
        'cf10812535b84e513d9dd75dc68f601915de10beb5326de42f2e1e88b1ae4054', ...
        '49f067005226bf21d61cb9c5f7d358c3a1cdead70c84ec93d637d38f068d7f43'});

results = cell(1, numel(routes));
for q = 1:numel(routes)
    source_file = fullfile(source_dir, routes(q).file);
    actual_sha = file_sha256(source_file);
    if ~strcmp(actual_sha, routes(q).sha256)
        error('ACI:FBCIROracleHash', ...
            ['Refusing to execute changed FBCIR source %s. Expected SHA-256 %s, ' ...
             'found %s. Audit the source before refreshing this harness.'], ...
             routes(q).file, routes(q).sha256, actual_sha);
    end
    results{q} = run_route(source_dir, source_file, output_dir, ...
                           routes(q).id, n_steps);
end

% The authors deliberately preserve the x1,x2,y1,y2 random-draw order in all
% three scripts. Check that the partitioned scripts really did generate the
% same physical path before shipping one shared signal file.
joint_signal = results{1}.signal;
for q = 2:numel(results)
    path_gap = max(abs(joint_signal(:) - results{q}.signal(:)));
    if path_gap > 5e-12
        error('ACI:FBCIROraclePath', ...
            'Partition %s did not reproduce the joint physical path (max gap %.17g).', ...
            routes(q).id, path_gap);
    end
end
signal_headers = {'t', 'x1', 'x2', 'y1', 'y2'};
write_numeric_csv(fullfile(output_dir, 'fbcir_multiscale_signal.csv'), ...
                  signal_headers, joint_signal);

record_path = fullfile(output_dir, 'fbcir_multiscale_generation_record.txt');
fid = fopen(record_path, 'w');
if fid < 0
    error('ACI:FBCIROracleWrite', 'Could not write %s.', record_path);
end
cleanup_record = onCleanup(@() fclose(fid));
fprintf(fid, 'fixture_kind: authors-source execution (not author-provided output)\n');
fprintf(fid, 'matlab_version: %s\n', version);
fprintf(fid, 'matlab_release: %s\n', version('-release'));
fprintf(fid, 'platform: %s\n', computer);
fprintf(fid, 'n_steps: %d\n', n_steps);
fprintf(fid, 'dt: %.17g\n', results{1}.dt);
fprintf(fid, 'seed_from_source: 444\n');
for q = 1:numel(routes)
    fprintf(fid, '%s_file: %s\n', routes(q).id, routes(q).file);
    fprintf(fid, '%s_sha256: %s\n', routes(q).id, routes(q).sha256);
end
clear cleanup_record

fprintf('Wrote FBCIR multiscale authors-source execution fixtures to %s\n', output_dir);
end


function result = run_route(source_dir, source_file, output_dir, route, n_steps)
text = fileread(source_file);
text = replace_once(text, 'N = 20000;', sprintf('N = %d;', n_steps), source_file);
text = replace_once(text, 'time_start_plot = 50;', 'time_start_plot = 0;', source_file);
text = replace_once(text, ...
    'calculate_and_plot_normalised_forward_CIR_metric = true;', ...
    'calculate_and_plot_normalised_forward_CIR_metric = false;', source_file);
text = replace_once(text, ...
    'calculate_and_plot_normalised_backward_CIR_metric = true;', ...
    'calculate_and_plot_normalised_backward_CIR_metric = false;', source_file);
text = strrep(text, 'figure(''WindowState'', ''maximized'');', ...
                    'figure(''Visible'', ''off'');');
text = replace_once(text, ...
    'clear online_fixed_mean online_fixed_cov update_matrices_fixed', ...
    '% fixture harness retains online_fixed_* through export', source_file);
text = replace_once(text, ...
    'clear online_adapt_mean online_adapt_cov update_matrices_adapt', ...
    '% fixture harness retains online_adapt_* through export', source_file);

temp_dir = tempname;
mkdir(temp_dir);
cleanup_temp = onCleanup(@() remove_temp_dir(temp_dir));
temp_script = fullfile(temp_dir, sprintf('fbcir_%s_fixture_source.m', route));
fid = fopen(temp_script, 'w');
if fid < 0
    error('ACI:FBCIROracleWrite', 'Could not write temporary source script.');
end
fprintf(fid, '%s', text);
fclose(fid);

old_visibility = get(groot, 'defaultFigureVisible');
visibility_cleanup = onCleanup(@() set(groot, 'defaultFigureVisible', old_visibility));
set(groot, 'defaultFigureVisible', 'off');
addpath(source_dir);
path_cleanup = onCleanup(@() rmpath(source_dir));

% RUN executes the hash-verified, minimally patched authors' script in this
% function workspace. The variables below are produced by that source.
run(temp_script);
close all force;

N1 = N + 1;
if N ~= n_steps || length(time) ~= N1
    error('ACI:FBCIROracleRun', 'The scale patch did not take effect for %s.', route);
end
sample_idx = unique(round(linspace(1, N1, min(N1, 61))));

if strcmp(route, 'joint')
    physical_signal = [time(:), x(1, :)', x(2, :)', y(1, :)', y(2, :)'];
else
    if strcmp(route, 'y1')
        physical_signal = [time(:), x(1, :)', x(2, :)', y(:), x(3, :)'];
    else
        physical_signal = [time(:), x(1, :)', x(2, :)', x(3, :)', y(:)];
    end
end

posterior_data = posterior_rows(sample_idx, time, filter_mean, filter_cov, ...
                                smoother_mean, smoother_cov, ACI_metric, ...
                                plot_idx, state_dim);
posterior_headers = posterior_header(state_dim);
write_numeric_csv(fullfile(output_dir, sprintf( ...
    'fbcir_multiscale_%s_posterior.csv', route)), ...
    posterior_headers, posterior_data);

pairs = online_pairs(N1);
fixed_data = online_rows(pairs, online_fixed_mean, online_fixed_cov, state_dim);
write_numeric_csv(fullfile(output_dir, sprintf( ...
    'fbcir_multiscale_%s_online_fixed.csv', route)), ...
    online_header(state_dim), fixed_data);

adapt_data = online_rows(pairs, online_adapt_mean, online_adapt_cov, state_dim);
write_numeric_csv(fullfile(output_dir, sprintf( ...
    'fbcir_multiscale_%s_online_adaptive.csv', route)), ...
    online_header(state_dim), adapt_data);

thresholds = [1e-1, 1e-2, 1e-3, 1e-4];
threshold_rows = zeros(size(thresholds));
for z = 1:numel(thresholds)
    [~, threshold_rows(z)] = min(abs(CIR_epsilon_values - thresholds(z)));
end
forward_data = [ ...
    sample_idx(:), time(sample_idx(:))', ...
    max_forw_RE_metric(sample_idx(:))', ...
    forw_approx_objective_CIR(sample_idx(:))', ...
    forw_defn_objective_CIR(sample_idx(:))', ...
    forw_subjective_CIR(threshold_rows, sample_idx(:))'];
range_headers = {'index', 't', 'max_metric', 'objective_l1_linf', ...
                 'objective_definition', 'subjective_1e_1', ...
                 'subjective_1e_2', 'subjective_1e_3', 'subjective_1e_4'};
write_numeric_csv(fullfile(output_dir, sprintf( ...
    'fbcir_multiscale_%s_forward_cir.csv', route)), ...
    range_headers, forward_data);

backward_data = [ ...
    sample_idx(:), time(sample_idx(:))', ...
    back_max_RE_metric(sample_idx(:))', ...
    back_approx_objective_CIR(sample_idx(:))', ...
    back_defn_objective_CIR(sample_idx(:))', ...
    back_subjective_CIR(threshold_rows, sample_idx(:))'];
write_numeric_csv(fullfile(output_dir, sprintf( ...
    'fbcir_multiscale_%s_backward_cir.csv', route)), ...
    range_headers, backward_data);

result = struct('signal', physical_signal, 'dt', dt);
clear path_cleanup visibility_cleanup cleanup_temp
end


function data = posterior_rows(idx, time, fm, fc, sm, sc, aci_metric, plot_idx, l)
N1 = length(time);
fm = reshape(fm, l, N1)';
sm = reshape(sm, l, N1)';
fc_flat = covariance_rows(fc, l, N1);
sc_flat = covariance_rows(sc, l, N1);
aci_full = nan(N1, 1);
aci_full(plot_idx) = aci_metric(:);
data = [idx(:), time(idx(:))', fm(idx, :), fc_flat(idx, :), ...
        sm(idx, :), sc_flat(idx, :), aci_full(idx)];
end


function headers = posterior_header(l)
headers = {'index', 't'};
for a = 1:l
    headers{end + 1} = sprintf('filter_mean_%d', a); %#ok<AGROW>
end
headers = [headers, covariance_header('filter_cov', l)];
for a = 1:l
    headers{end + 1} = sprintf('smoother_mean_%d', a); %#ok<AGROW>
end
headers = [headers, covariance_header('smoother_cov', l), {'ACI_metric'}];
end


function pairs = online_pairs(N1)
j_values = unique([1, 2, round(N1 / 5), round(N1 / 2), N1 - 1, N1]);
pairs = zeros(0, 2);
for j = j_values
    n_values = unique([j, min(N1, j + 1), round((j + N1) / 2), N1]);
    n_values = n_values(n_values >= j & n_values <= N1);
    pairs = [pairs; [repmat(j, numel(n_values), 1), n_values(:)]]; %#ok<AGROW>
end
pairs = unique(pairs, 'rows', 'stable');
end


function data = online_rows(pairs, means, covariances, l)
data = zeros(size(pairs, 1), 2 + l + l * l);
data(:, 1:2) = pairs;
for q = 1:size(pairs, 1)
    j = pairs(q, 1);
    n = pairs(q, 2);
    mu = reshape(means{n}(:, j), 1, l);
    R = reshape(covariances{n}(:, :, j), l, l);
    data(q, 3:(2 + l)) = mu;
    data(q, (3 + l):end) = reshape(R, 1, l * l);
end
end


function headers = online_header(l)
headers = {'j', 'n'};
for a = 1:l
    headers{end + 1} = sprintf('online_mean_%d', a); %#ok<AGROW>
end
headers = [headers, covariance_header('online_cov', l)];
end


function flat = covariance_rows(covariances, l, N1)
flat = zeros(N1, l * l);
if l == 1
    flat(:, 1) = reshape(covariances, N1, 1);
    return
end
for q = 1:N1
    R = reshape(covariances(:, :, q), l, l);
    flat(q, :) = reshape(R, 1, l * l);
end
end


function headers = covariance_header(prefix, l)
headers = cell(1, l * l);
z = 0;
% MATLAB reshape order is column-major: (1,1),(2,1),...,(l,l).
for b = 1:l
    for a = 1:l
        z = z + 1;
        headers{z} = sprintf('%s_%d_%d', prefix, a, b);
    end
end
end


function text = replace_once(text, old, new, source_file)
hits = strfind(text, old);
if numel(hits) ~= 1
    error('ACI:FBCIROraclePatch', ...
        'Expected one occurrence of "%s" in %s, found %d.', ...
        old, source_file, numel(hits));
end
text = strrep(text, old, new);
end


function write_numeric_csv(path, headers, values)
if size(values, 2) ~= numel(headers)
    error('ACI:FBCIROracleWrite', 'Header/data width mismatch for %s.', path);
end
fid = fopen(path, 'w');
if fid < 0
    error('ACI:FBCIROracleWrite', 'Could not write %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', headers{1});
for q = 2:numel(headers)
    fprintf(fid, ',%s', headers{q});
end
fprintf(fid, '\n');
for i = 1:size(values, 1)
    fprintf(fid, '%.17g', values(i, 1));
    for q = 2:size(values, 2)
        fprintf(fid, ',%.17g', values(i, q));
    end
    fprintf(fid, '\n');
end
clear cleanup
end


function hash = file_sha256(path)
fid = fopen(path, 'r');
if fid < 0
    error('ACI:FBCIROracleSource', 'Could not read %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
digest = typecast(engine.digest(), 'uint8');
hash = lower(reshape(dec2hex(digest, 2)', 1, []));
clear cleanup
end


function remove_temp_dir(path)
if isfolder(path)
    rmdir(path, 's');
end
end
