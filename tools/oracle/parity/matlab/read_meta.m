function meta = read_meta(path)
%READ_META  Parse a dataset's meta.dcf into a struct.
%
%   Debian-control format: `Field: value` per line, `#` comments ignored,
%   continuation lines (leading whitespace) appended. Numeric-looking values
%   become doubles, everything else stays a string. Twenty lines rather than a
%   toolbox dependency, and it reads the same file R does.

meta = struct();
fid = fopen(path, 'r');
if fid < 0
    error('parity:noMeta', 'cannot open %s', path);
end
cleaner = onCleanup(@() fclose(fid));

field = '';
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    if isempty(strtrim(line)) || startsWith(strtrim(line), '#')
        continue
    end
    if ~isempty(regexp(line, '^\s', 'once')) && ~isempty(field)
        meta.(field) = [meta.(field) ' ' strtrim(line)];
        continue
    end
    parts = regexp(line, '^([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(.*)$', 'tokens', 'once');
    if isempty(parts)
        continue
    end
    field = parts{1};
    meta.(field) = strtrim(parts{2});
end

names = fieldnames(meta);
for i = 1:numel(names)
    value = str2double(meta.(names{i}));
    if ~isnan(value)
        meta.(names{i}) = value;
    end
end

end
