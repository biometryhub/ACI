function probe_conditional(workspaces, profiles, out_path)
%PROBE_CONDITIONAL  Measure whether a captured ENSO run assimilated one channel.
%
%   Finding F7 asked which of the five ENSO scripts ship with conditional ACI
%   enabled, and answered it by reading the sources. Reading is how the
%   question was got wrong twice. This measures it instead, from the workspaces
%   the scripts themselves produced.
%
%   The discriminator is the observation-noise inverse the filter actually
%   used. Two structures are possible and they are not close:
%
%     unconditional -- the per-step pseudoinverse is live, so S_xoS_x_inv
%                      varies from slice to slice and is generally dense;
%     conditional   -- the per-step pseudoinverse is commented out and one
%                      diagonal entry is overwritten across all slices, so
%                      every slice from the second onward is the SAME matrix
%                      with exactly one non-zero entry.
%
%   Slice one is reported separately, because it is the one slice the
%   conditional scripts leave partly dense: the pseudoinverse fills it before
%   the single entry is overwritten, and the filter consumes it at j = 2.

rows = cell(numel(workspaces), 1);
for k = 1:numel(workspaces)
    w = load(workspaces{k}, 'S_xoS_x_inv');
    A = w.S_xoS_x_inv;
    n = size(A, 3);
    d = size(A, 1);

    first = A(:, :, 1);
    second = A(:, :, 2);

    % Are all slices from the second onward identical? A live per-step
    % pseudoinverse makes them differ.
    tail_constant = true;
    tail_spread = 0;
    for j = 3:n
        difference = max(abs(A(:, :, j) - second), [], 'all');
        tail_spread = max(tail_spread, difference);
        if difference > 0
            tail_constant = false;
        end
    end

    nz_first = nnz(first);
    nz_second = nnz(second);

    % The entries slice one carries that the later slices do not -- the
    % pseudoinverse residue a conditional run assimilates once and never again.
    off = first;
    off(second ~= 0) = 0;
    first_extra = max(abs(off), [], 'all');

    if tail_constant && nz_second == 1
        verdict = 'conditional';
    elseif ~tail_constant
        verdict = 'unconditional';
    else
        verdict = 'indeterminate';
    end

    rows{k} = {profiles{k}, d, n, nz_first, nz_second, ...
               double(tail_constant), ...
               tail_spread, first_extra, second(1, 1), verdict};
    fprintf(['%-14s %dx%d x%d : slice1 nnz %d, slice2 nnz %d, tail ' ...
             'constant %d (spread %.3g), slice1 extra %.6g, (1,1) ' ...
             '%.10g -> %s\n'], ...
            profiles{k}, d, d, n, nz_first, nz_second, tail_constant, ...
            tail_spread, first_extra, second(1, 1), verdict);
end

T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'profile', 'dim', 'slices', 'nnz_slice1', 'nnz_slice2', ...
     'tail_constant', 'tail_spread', 'slice1_extra', 'entry_11', 'verdict'});
writetable(T, out_path);
fprintf('wrote %s\n', out_path);

end
