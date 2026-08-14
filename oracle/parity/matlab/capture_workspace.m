function capture_workspace(script_path, reference_dir, out_path)
%CAPTURE_WORKSPACE  Run a reduced reference script and save its workspace.
%
%   The gold side of gate G1. The script is executed unaltered except for the
%   problem-size substitutions applied by tools/reduce.R, whose diff against
%   the reference is three lines and is reproduced in the parity report. Every
%   variable the script leaves behind is saved, so an extracted function can
%   later be checked against the script's own values rather than against a
%   transcription of them.
%
%   Figures are suppressed rather than skipped: the plotting sections are part
%   of the script and are allowed to run, but a plotting failure must not cost
%   us the computed workspace, so the run is wrapped and the save happens
%   either way. Whether it errored is recorded in `capture_error`.

set(0, 'DefaultFigureVisible', 'off');
addpath(reference_dir);                    % progress_bar, simps, legendUnq
addpath(fileparts(script_path));

capture_error = '';
capture_started = tic;
try
    run(script_path);
catch capture_ME
    capture_error = sprintf('%s: %s', capture_ME.identifier, capture_ME.message);
    fprintf(2, 'script raised: %s\n', capture_error);
end
capture_seconds = toc(capture_started);

close all force;
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end
save(out_path, '-v7.3');
fprintf('captured %s (%.1f s)\n', out_path, capture_seconds);

end
