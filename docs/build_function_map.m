function build_function_map(root, outdir)
%BUILD_FUNCTION_MAP Rebuilds function_map.html from the toolbox code.
%
%   BUILD_FUNCTION_MAP
%   BUILD_FUNCTION_MAP(root)
%   BUILD_FUNCTION_MAP(root, outdir)
%
%   Scans every .m file in bvartools/, cmintools/ and examples/, finds which
%   toolbox functions each file calls, works out which functions are used
%   by the examples, reads the help text of each function and writes the
%   page docs/function_map.html (dependency graph and function register).
%   Run it after adding, removing or renaming functions, or after editing
%   their help, and commit the new page.
%
%   Inputs
%     root     folder of the toolbox (default: the parent of this file's folder)
%     outdir   folder where function_map.html is written, or a cell array of
%              folders to write it to several places (default: root/docs).
%              The folders must exist. Example, toolbox and personal site:
%                build_function_map('C:\Users\Filippo\Documents\BVAR_', ...
%                  {'C:\Users\Filippo\Documents\BVAR_\docs', ...
%                   'C:\Users\Filippo\Documents\naffe15.github.io\_pages'})
%
%   How a call is detected
%     - comments and strings are removed from each line before matching;
%     - a name counts as a call when it is the name of a toolbox function and
%       the calling file does not use it as a variable or define a local
%       function with that name;
%     - optimiser objectives passed by name ('bvar_opt_hyperpara' etc.) count
%       as calls: any string equal to a toolbox function containing '_opt_';
%     - calls built with eval or feval of a variable are not detected.
%
%   Status of each function
%     entry      called by at least one example script
%     internal   reached by following the calls out from the entry points
%     dead       called only by functions that are never reached
%     orphan     never called
%
%   The "Guide" column counts calls like name( in HitchhikerGuide_.pdf. It is
%   filled only if the pdftotext program (Xpdf/Poppler) is installed;
%   otherwise it is left empty.
%
%   Files
%     docs/function_map_template.html   page layout; the data are inserted
%                                       at /*__DATA__*/ and /*__META__*/
%     function_map.html                 the page written by this function,
%                                       in each folder of outdir
%
%   Works in MATLAB (R2016b or later) and GNU Octave.

%   Filippo Ferroni and Fabio Canova. https://github.com/naffe15/BVAR_

if nargin < 1 || isempty(root)
    root = fileparts(fileparts(mfilename('fullpath')));
end
docs = fullfile(root, 'docs');                % the template is read from here
if nargin < 2 || isempty(outdir)
    outdir = {docs};
elseif ischar(outdir)
    outdir = {outdir};
end
for k = 1 : numel(outdir)
    if exist(outdir{k}, 'dir') ~= 7
        error('build_function_map: output folder not found: %s', outdir{k});
    end
end

% ------------------------------------------------------------------ files
groups = {'bvartools', 'cmintools', 'examples'};
F = struct('path', {}, 'group', {}, 'name', {}, 'rel', {});
for g = 1 : numel(groups)
    list = list_m_files(fullfile(root, groups{g}));
    for k = 1 : numel(list)
        [~, nm] = fileparts(list{k});
        F(end+1).path = list{k};                                   %#ok<AGROW>
        F(end).group  = groups{g};
        F(end).name   = nm;
        F(end).rel    = strrep(list{k}(numel(root)+2 : end), '\', '/');
    end
end
istool = ~strcmp({F.group}, 'examples');
tool   = containers.Map({F(istool).name}, {F(istool).group});
fprintf('build_function_map: %d toolbox functions, %d example scripts\n', ...
    sum(istool), sum(~istool));

% ------------------------------------------------------------ scan files
S = cell(numel(F), 1);
for k = 1 : numel(F)
    S{k} = scan_file(F(k).path);
end

% ------------------------------------------------------------ call graph
calls   = containers.Map('KeyType', 'char', 'ValueType', 'any');  % caller -> callees
callers = containers.Map('KeyType', 'char', 'ValueType', 'any');  % callee -> callers
for k = 1 : numel(F)
    if istool(k)
        src = F(k).name;
    else
        src = ['EX:' F(k).rel];
    end
    s   = S{k};
    out = {};
    for r = 1 : numel(s.refs)
        x = s.refs{r};
        if isKey(tool, x) && ~strcmp(x, F(k).name) && ~any(strcmp(x, s.defs)) ...
                && ~any(strcmp(x, s.vars))
            out{end+1} = x;                                         %#ok<AGROW>
        end
    end
    for r = 1 : numel(s.strs)
        x = s.strs{r};
        if isKey(tool, x) && ~isempty(strfind(x, '_opt_')) && ~strcmp(x, F(k).name)
            out{end+1} = x;                                         %#ok<AGROW>
        end
    end
    out = unique(out);
    calls(src) = out;
    for r = 1 : numel(out)
        if isKey(callers, out{r})
            callers(out{r}) = [callers(out{r}) {src}];
        else
            callers(out{r}) = {src};
        end
    end
end

% entry points: functions called by the examples; then follow the calls
roots = {};
for k = find(~istool)
    roots = [roots calls(['EX:' F(k).rel])];                        %#ok<AGROW>
end
roots = unique(roots);
reach = {};
stack = roots;
while ~isempty(stack)
    x = stack{end};
    stack(end) = [];
    if any(strcmp(x, reach))
        continue
    end
    reach{end+1} = x;                                               %#ok<AGROW>
    if isKey(calls, x)
        stack = [stack calls(x)];                                   %#ok<AGROW>
    end
end

% ------------------------------------------------- Guide mentions (optional)
guide = '';
pdf   = fullfile(root, 'HitchhikerGuide_.pdf');
if exist(pdf, 'file') == 2
    [st, ~] = system('pdftotext -v 2>&1');
    if st == 0
        txt = [tempname '.txt'];
        [st, ~] = system(['pdftotext -layout "' pdf '" "' txt '" 2>&1']);
        if st == 0
            guide = fileread(txt);
            delete(txt);
        end
    end
end
if isempty(guide)
    fprintf('build_function_map: pdftotext not found, the Guide column is left empty\n');
end

% ---------------------------------------------------------------- rows
% [name, folder, status, examples, callers, guide, lines, '', help, signature,
%  local functions, description, has references, help lines]
idx   = find(istool);
keys_ = cell(numel(idx), 1);
for j = 1 : numel(idx)
    keys_{j} = sprintf('%d|%s', ~strcmp(F(idx(j)).group, 'bvartools'), lower(F(idx(j)).name));
end
[~, ord] = sort(keys_);
idx  = idx(ord);
rows = cell(1, numel(idx));
nst  = struct('E', 0, 'I', 0, 'D', 0, 'O', 0);
for j = 1 : numel(idx)
    k  = idx(j);
    n  = F(k).name;
    s  = S{k};
    cl = {};
    if isKey(callers, n), cl = unique(callers(n)); end
    isex = strncmp(cl, 'EX:', 3);
    ex   = cellfun(@example_label, cl(isex), 'UniformOutput', false);
    tc   = cl(~isex);
    if any(strcmp(n, roots))
        st = 'E';
    elseif any(strcmp(n, reach))
        st = 'I';
    elseif ~isempty(tc)
        st = 'D';
    else
        st = 'O';
    end
    nst.(st) = nst.(st) + 1;
    % help: description from the first line, toolbox layout or third party
    ht   = strjoin(s.help, char(10));
    desc = '';
    if ~isempty(s.help)
        h1  = regexprep(s.help{1}, '^%\s*', '');
        tok = regexp(h1, '^([A-Za-z_0-9]+):?\s+(.*)$', 'tokens', 'once');
        if ~isempty(tok) && strcmpi(tok{1}, n)
            desc = tok{2};
        else
            desc = h1;
        end
    end
    if numel(desc) > 140, desc = desc(1:140); end
    third = ~isempty(regexpi(ht, 'original help below', 'once')) || ...
        any(strcmp(n, {'autocor', 'distinguishable_colors', 'colorspace'}));
    if third, hl = 't'; else, hl = 'f'; end
    hasref = ~isempty(regexpi(ht, 'Reference|Copyright|Original file|Adapted from|Based on|code by', 'once'));
    gd = 0;
    if ~isempty(guide)
        gd = numel(regexp(guide, ['(?<![\w.])' regexptranslate('escape', n) '(\s*\(|\.m(?!\w))']));
    end
    rows{j} = {n, F(k).group(1), st, ex, tc, gd, s.nlines, '', hl, s.sig, ...
        s.defs(2:end), desc, double(hasref), numel(s.help)};
end

% ---------------------------------------------------------------- write
commit = '';
[st, out] = system(['git -C "' root '" rev-parse --short HEAD 2>&1']);
if st == 0, commit = strtrim(out); end
meta = {datestr(now, 'dd mmm yyyy'), commit, sum(~istool), ~isempty(guide)};
page = fileread(fullfile(docs, 'function_map_template.html'));
if isempty(strfind(page, '/*__DATA__*/')) || isempty(strfind(page, '/*__META__*/'))
    error('build_function_map: the template has no /*__DATA__*/ or /*__META__*/ marker.');
end
page = strrep(page, '/*__DATA__*/', strrep(tojson(rows), '</', '<\/'));
page = strrep(page, '/*__META__*/', strrep(tojson(meta), '</', '<\/'));
fprintf('build_function_map: %d entry, %d internal, %d dead, %d orphan\n', ...
    nst.E, nst.I, nst.D, nst.O);
for k = 1 : numel(outdir)
    out = fullfile(outdir{k}, 'function_map.html');
    fid = fopen(out, 'w');
    if fid < 0
        error('build_function_map: cannot write %s', out);
    end
    fwrite(fid, page);
    fclose(fid);
    fprintf('build_function_map: written %s\n', out);
end
end


% =========================================================================
function list = list_m_files(folder)
% All .m files in folder and its subfolders.
list = {};
d    = dir(folder);
for k = 1 : numel(d)
    if d(k).isdir
        if ~any(strcmp(d(k).name, {'.', '..'}))
            list = [list list_m_files(fullfile(folder, d(k).name))];   %#ok<AGROW>
        end
    elseif numel(d(k).name) > 2 && strcmp(d(k).name(end-1:end), '.m')
        list{end+1} = fullfile(folder, d(k).name);                    %#ok<AGROW>
    end
end
end


% =========================================================================
function s = scan_file(path)
% Definitions, variables, identifiers and strings of one file, plus its help.
txt   = fileread(path);
txt(double(txt) > 126) = '?';                % non-ASCII characters (only found in
                                             % comments) would stop Octave's regexp
lines = regexp(txt, '\r?\n', 'split');
s.nlines = numel(lines);
s.defs = {};  s.vars = {};  s.refs = {};  s.strs = {};
inblock = false;
skip    = 0;                                 % last line already read as a continuation
for k = 1 : numel(lines)
    if k <= skip, continue, end
    raw = lines{k};
    t   = strtrim(raw);
    if inblock                               % inside a %{ ... %} block
        if strcmp(t, '%}'), inblock = false; end
        continue
    end
    if strcmp(t, '%{'), inblock = true; continue, end
    [code, strs] = strip_line(raw);
    if isempty(strtrim(code)), continue, end
    % function definition: its name and arguments
    if ~isempty(regexp(code, '^\s*function[\s\[]', 'once'))
        while ~isempty(regexp(code, '\.\.\.\s*$', 'once')) && k < numel(lines)
            k    = k + 1;                    % function line continued with ...
            code = [regexprep(code, '\.\.\.\s*$', ' ') strip_line(lines{k})];
            skip = k;
        end
        [nm, args] = parse_function_line(code);
        s.defs{end+1} = nm;                                          %#ok<AGROW>
        s.vars = [s.vars args];                                      %#ok<AGROW>
        continue
    end
    % variables: [a,b] = ..., a = ..., a(i).f = ..., for a = ...
    tok = regexp(code, '^\s*\[([^\]=]*)\]\s*=(?!=)', 'tokens', 'once');
    if ~isempty(tok)
        parts = regexp(tok{1}, '[\s,]+', 'split');
        for p = 1 : numel(parts)
            v = regexp(parts{p}, '^[A-Za-z]\w*', 'match', 'once');
            if ~isempty(v), s.vars{end+1} = v; end                   %#ok<AGROW>
        end
    end
    tok = regexp(code, '^\s*([A-Za-z]\w*)\s*(\(.*?\)|\{.*?\})?\s*(\.\w+)*\s*=(?!=)', 'tokens', 'once');
    if ~isempty(tok), s.vars{end+1} = tok{1}; end                    %#ok<AGROW>
    tok = regexp(code, '^\s*for\s+([A-Za-z]\w*)\s*=', 'tokens', 'once');
    if ~isempty(tok), s.vars{end+1} = tok{1}; end                    %#ok<AGROW>
    % identifiers and identifier-like strings
    s.refs = [s.refs regexp(code, '(?<![\w.])[A-Za-z]\w*', 'match')];  %#ok<AGROW>
    for p = 1 : numel(strs)
        if ~isempty(regexp(strtrim(strs{p}), '^[A-Za-z]\w*$', 'once'))
            s.strs{end+1} = strtrim(strs{p});                        %#ok<AGROW>
        end
    end
end
s.vars = unique(s.vars);
s.refs = unique(s.refs);  s.strs = unique(s.strs);
% help as MATLAB reads it: the first unbroken comment block after the
% function line (which may continue over several lines with ...)
s.help = {};
s.sig  = '(script)';
f = find(~cellfun(@isempty, regexp(lines, '^\s*function[\s\[]', 'once')), 1);
if ~isempty(f)
    e = f;
    while e < numel(lines) && ~isempty(regexp(lines{e}, '\.\.\.\s*$', 'once'))
        e = e + 1;
    end
    sig = strjoin(cellfun(@strtrim, lines(f:e), 'UniformOutput', false), ' ');
    s.sig = regexprep(sig, '\s*\.\.\.\s*', ' ');
    i = e + 1;
    while i <= numel(lines) && isempty(strtrim(lines{i}))
        i = i + 1;
    end
    while i <= numel(lines) && ~isempty(regexp(lines{i}, '^\s*%', 'once'))
        s.help{end+1} = strtrim(lines{i});                           %#ok<AGROW>
        i = i + 1;
    end
end
end


% =========================================================================
function [nm, args] = parse_function_line(code)
% Name, outputs and inputs of 'function [a,b] = name(c,d)'. Done in steps
% because MATLAB and Octave return optional regexp groups differently.
rest = regexprep(code, '^\s*function\s*', '');
p    = strfind(rest, '(');
q    = strfind(rest, '=');
outs = '';
if ~isempty(q) && (isempty(p) || q(1) < p(1))
    outs = rest(1 : q(1)-1);
    rest = rest(q(1)+1 : end);
end
nm   = regexp(rest, '^\s*(\w+)', 'tokens', 'once');
nm   = nm{1};
ins  = regexp(rest, '\((.*)\)', 'tokens', 'once');
if isempty(ins), ins = {''}; end
args = regexp([outs ',' ins{1}], '[A-Za-z]\w*', 'match');
end


% =========================================================================
function [code, strs] = strip_line(raw)
% Removes the comment and the string literals from one line of code. Each
% string is replaced by a blank and returned in strs. A quote is a transpose
% when it follows a letter, digit, _, ), ], }, . or another transpose.
strs = {};
if isempty(regexp(raw, '[%''"]', 'once'))
    code = raw;
    return
end
idc = ['A':'Z' 'a':'z' '0':'9' '_)]}.'''];   % characters a transpose can follow
out = blanks(numel(raw));
n   = 0;
i   = 1;
L   = numel(raw);
while i <= L
    c = raw(i);
    if c == '%'
        break
    end
    istr = c == '"' || (c == '''' && ~(n > 0 && any(out(n) == idc)));
    if istr
        j   = i + 1;
        buf = '';
        while j <= L
            if raw(j) == c
                if j < L && raw(j+1) == c        % doubled quote inside the string
                    buf(end+1) = c;              %#ok<AGROW>
                    j = j + 2;
                    continue
                end
                break
            end
            buf(end+1) = raw(j);                 %#ok<AGROW>
            j = j + 1;
        end
        strs{end+1} = buf;                       %#ok<AGROW>
        n = n + 1;
        out(n) = ' ';
        i = j + 1;
        continue
    end
    n = n + 1;
    out(n) = c;
    i = i + 1;
end
code = out(1:n);
end


% =========================================================================
function lab = example_label(src)
% 'EX:examples/BVAR tutorial/example_3_irf.m' -> 'bvar/example_3_irf'
parts = regexp(src(4:end), '/', 'split');
folder = parts{end-1};
switch folder
    case 'BVAR tutorial Octave',        folder = 'oct';
    case 'BVAR tutorial',               folder = 'bvar';
    case 'Trend-Cycle-Dating tutorial', folder = 'tcd';
end
[~, nm] = fileparts(parts{end});
lab = [folder '/' nm];
end


% =========================================================================
function j = tojson(x)
% Minimal JSON writer for cell arrays (lists), char (strings), numbers and
% logicals. Non-ASCII characters are written as \uXXXX, so the page is
% plain ASCII whatever the encoding of the source files.
if iscell(x)
    parts = cellfun(@tojson, x, 'UniformOutput', false);
    j = ['[' strjoin(parts, ',') ']'];
elseif ischar(x)
    j = '"';
    for c = double(x(:)')
        if c == 34
            j = [j '\"'];                                            %#ok<AGROW>
        elseif c == 92
            j = [j '\\'];                                            %#ok<AGROW>
        elseif c < 32 || c > 126
            j = [j sprintf('\\u%04x', c)];                           %#ok<AGROW>
        else
            j = [j char(c)];                                         %#ok<AGROW>
        end
    end
    j = [j '"'];
elseif islogical(x) || isnumeric(x)
    if isscalar(x)
        j = num2str(double(x));
    else
        j = tojson(num2cell(double(x(:)')));
    end
else
    error('tojson: unsupported type %s', class(x));
end
end
