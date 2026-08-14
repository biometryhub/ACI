function failures = check_profile(workspace_path, conditions)
%CHECK_PROFILE  Assert a captured workspace is not degenerate.
%
%   A reduced problem size can quietly hollow out the very quantity it was
%   meant to grade. The reference indexes its exact objective CIR as
%   `subjective_CIR(:, 1:end-lookahead_tolerance/dt)`; when the reporting
%   window is shorter than the lookahead, MATLAB evaluates that to an EMPTY
%   range rather than raising, and the quantity silently becomes 1x0. A pairing
%   graded against it would then compare nothing to nothing and pass.
%
%   CONDITIONS is a cell array of expression strings declared in
%   manifest/knobs.dcf under `Requires`. Each is evaluated against the loaded
%   workspace and must return a scalar true.
%
%   Returns the failing conditions, and prints each verdict.

w = load(workspace_path);
names = fieldnames(w);
for i = 1:numel(names)
    eval([names{i} ' = w.(names{i});']);       %#ok<EVLDIR> -- declared checks
end

failures = {};
for i = 1:numel(conditions)
    expression = strtrim(conditions{i});
    if isempty(expression)
        continue
    end
    try
        verdict = eval(expression);
        passed = islogical(verdict) && isscalar(verdict) && verdict;
    catch ME
        passed = false;
        fprintf(2, '  ERROR  %s  (%s)\n', expression, ME.message);
    end
    if passed
        fprintf('  ok     %s\n', expression);
    else
        fprintf(2, '  FAIL   %s\n', expression);
        failures{end+1} = expression;             %#ok<AGROW>
    end
end

if isempty(failures)
    fprintf('profile check passed (%d conditions)\n', numel(conditions));
else
    error('parity:degenerateProfile', ...
          '%d of %d profile conditions failed for %s', ...
          numel(failures), numel(conditions), workspace_path);
end

end
