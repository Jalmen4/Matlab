%% plot_results_ccs_script.m — data-driven loader (Transfer first, gdxdump fallback)
% Robustly builds arrays directly from symbol records (no reliance on set sizes).
% Final shapes:
%   qCO2 : [Q × N × N × T]
%   cCO2 : [C × N × T]
%   Seq  : [N × T]
% Saves: ./figures/costs.png, map.png, sankey.png

close all; clc;

%% ---- 0) Paths ------------------------------------------------------------
% If needed, make sure the GAMS Transfer API is on your path:
% addpath('C:\GAMS\48\api\matlab');   % <-- uncomment & adjust if required

thisFile = mfilename('fullpath'); if isempty(thisFile), thisFile = which('plot_results_ccs_script.m'); end
projectRoot = fileparts(fileparts(thisFile));
figDir = fullfile(projectRoot, 'figures'); if ~exist(figDir, 'dir'), mkdir(figDir); end

% === EDIT ONE OF THESE LINES TO POINT TO YOUR Results.gdx ==================
% gdxResults = fullfile(projectRoot, 'data', 'Results.gdx');             % SPEC default
gdxResults = fullfile(projectRoot, 'Programa nuevo 1.0', 'Results.gdx'); % your custom path
% ===========================================================================
assert(exist(gdxResults,'file')==2, 'Results.gdx not found at: %s', gdxResults);

%% ---- 1) Load symbol *tables* (Transfer if possible; else gdxdump) -------
use_gdxdump = true;
R_qCO2 = table(); R_cCO2 = table(); R_Seq = table();
R_TCC  = table(); R_TTC  = table(); R_TSC = table(); R_TC = table(); R_Z = table();
R_coord= table(); R_coastal = table(); R_mar = table(); R_inland = table();
R_aco  = table(); R_UTC = table(); R_UCC = table(); R_USC = table(); R_USCoff = table();

% Helper to normalize a records payload into a table with named columns
symbolToTable = @(sym, names) normalizeSymbolRecords(sym, names);

if exist('gams.transfer.Container','class')==8
    import gams.transfer.*
    try
        try m = Container.read(gdxResults); catch, m = Container(); m.read(gdxResults); end
        % Try to fetch the three core value symbols; if any fails → fallback to gdxdump
        ok = true;
        try R_qCO2 = symbolToTable(fetchSym(m,'qCO2_l'), {'q','n','nn','t','value'}); catch, ok = false; end
        try R_cCO2 = symbolToTable(fetchSym(m,'cCO2_l'), {'c','n','t','value'});      catch, ok = false; end
        try R_Seq  = symbolToTable(fetchSym(m,'Seq_l'),  {'n','t','value'});          catch, ok = false; end
        if ok
            % Optional/aux symbols (best-effort)
            R_TCC   = trySymTable(m,'TCC_l', {'t','value'});
            R_TTC   = trySymTable(m,'TTC_l', {'t','value'});
            R_TSC   = trySymTable(m,'TSC_l', {'t','value'});
            R_TC    = trySymTable(m,'TC_l',  {'t','value'});
            R_Z     = trySymTable(m,'Z_l',   {'value'});
            R_coord = trySymTable(m,'coord', {'n','p','value'});
            R_coastal = trySymTable(m,'coastal_n', {'n','value'});
            R_mar     = trySymTable(m,'mar_n',     {'n','value'});
            R_inland  = trySymTable(m,'inland_n',  {'n','value'});
            R_aco     = trySymTable(m,'aco',       {'n','nn','value'});
            R_UTC     = trySymTable(m,'UTC',       {'q','value'});
            R_UCC     = trySymTable(m,'UCC',       {'c','value'});
            R_USC     = trySymTable(m,'USC',       {'value'});
            R_USCoff  = trySymTable(m,'USC_offshore', {'n','value'});
            use_gdxdump = false;
        end
    catch
        use_gdxdump = true;
    end
end

if use_gdxdump
    % --- CLI fallback using gdxdump (auto-detect path) ---
    gdxResultsQ = ['"' gdxResults '"'];
    gdxdumpCmd = autodetect_gdxdump();   % throws if not found
    dumpCSV = @(sym, cols) readCsvDump(gdxdumpCmd, gdxResultsQ, sym, cols);

    R_qCO2  = dumpCSV('qCO2_l',  {'q','n','nn','t','value'});
    R_cCO2  = dumpCSV('cCO2_l',  {'c','n','t','value'});
    R_Seq   = dumpCSV('Seq_l',   {'n','t','value'});
    R_TCC   = dumpCSV('TCC_l',   {'t','value'});
    R_TTC   = dumpCSV('TTC_l',   {'t','value'});
    R_TSC   = dumpCSV('TSC_l',   {'t','value'});
    R_TC    = dumpCSV('TC_l',    {'t','value'});
    R_Z     = dumpCSV('Z_l',     {'value'});
    R_coord = dumpCSV('coord',   {'n','p','value'});
    R_coastal = dumpCSV('coastal_n', {'n','value'});
    R_mar     = dumpCSV('mar_n',     {'n','value'});
    R_inland  = dumpCSV('inland_n',  {'n','value'});
    R_aco     = dumpCSV('aco',       {'n','nn','value'});
    R_UTC     = dumpCSV('UTC',       {'q','value'});
    R_UCC     = dumpCSV('UCC',       {'c','value'});
    R_USC     = dumpCSV('USC',       {'value'});
    R_USCoff  = dumpCSV('USC_offshore', {'n','value'});
end

%% ---- 2) Data-driven domains (derive from records, not sets) -------------
% Get unique labels in the order they appear
uniq = @(x) unique(x, 'stable');
Qlab = uniq(toStrCol(R_qCO2.q));
Nlab = uniq(toStrCol(R_qCO2.n));
NNlab= uniq(toStrCol(R_qCO2.nn));
Tlab = uniq(toStrCol(R_qCO2.t));

% Sanity: n and nn domains should be the same set; if not, unify
if ~isequal(sort(Nlab), sort(NNlab))
    Nall = uniq([Nlab; NNlab]);
else
    Nall = Nlab;
end
Nlab = Nall; NNlab = Nall;

% Also derive C from cCO2
Clab = uniq(toStrCol(R_cCO2.c));

Q = numel(Qlab); N = numel(Nlab); C = numel(Clab); T = numel(Tlab);

%% ---- 3) Build dense arrays (robust to column order) ---------------------
qCO2 = tableToDense(R_qCO2, {'q','n','nn','t'}, {Qlab,Nlab,NNlab,Tlab}, [Q N N T]);
cCO2 = tableToDense(R_cCO2, {'c','n','t'},     {Clab,Nlab,Tlab},        [C N T]);
Seq  = tableToDense(R_Seq,  {'n','t'},         {Nlab,Tlab},             [N T]);

% Best-effort aux data (some files may not include all)
coord        = tableToDense(R_coord, {'n','p'}, {Nlab,{'Latitude','Longitude'}}, [N 2], true);
coastal_n    = tableToDense(R_coastal, {'n'}, {Nlab}, [N 1], true)>0.5;
mar_n        = tableToDense(R_mar,     {'n'}, {Nlab}, [N 1], true)>0.5;
inland_n     = tableToDense(R_inland,  {'n'}, {Nlab}, [N 1], true)>0.5;
aco          = tableToDense(R_aco,     {'n','nn'}, {Nlab,NNlab}, [N N], true);
UTC          = tableToDense(R_UTC,     {'q'}, {Qlab}, [Q 1], true);
UCC          = tableToDense(R_UCC,     {'c'}, {Clab}, [C 1], true);
USC          = scalarFromTable(R_USC);
USC_offshore = tableToDense(R_USCoff,  {'n'}, {Nlab}, [N 1], true);

% Time series
TCC = vecFromTable(R_TCC, {'t'}, {Tlab}, [T 1]);
TTC = vecFromTable(R_TTC, {'t'}, {Tlab}, [T 1]);
TSC = vecFromTable(R_TSC, {'t'}, {Tlab}, [T 1]);
TCv = vecFromTable(R_TC,  {'t'}, {Tlab}, [T 1]);  % 'TC' reserved word risk
Zl  = scalarFromTable(R_Z);

fprintf('[Loader] q>0=%d  c>0=%d  Seq>0=%d  |  T=%d N=%d Q=%d C=%d\n', ...
    nnz(qCO2>0), nnz(cCO2>0), nnz(Seq>0), T, N, Q, C);

%% ---- 4) Derived tensors --------------------------------------------------
qCO2_total = squeeze(sum(qCO2,1));  % N×N×T
lat = coord(:,1); lon = coord(:,2);
cls = zeros(N,1); cls(inland_n)=1; cls(coastal_n)=2; cls(mar_n)=3;

%% ---- 5) Costs figure -----------------------------------------------------
f1 = figure('Color','w','Name','CCS Costs'); ax1 = axes(f1); hold(ax1,'on');
tt = (1:T)';
bar(ax1, tt, [TCC(:) TTC(:) TSC(:)], 'stacked');
plot(ax1, tt, TCv(:), 'k.-', 'LineWidth', 1.2);
grid(ax1,'on'); xlabel(ax1,'Year (t index)'); ylabel(ax1,'Cost (€/yr)');
legend({'Capture','Transport','Sequestration','Total'}, 'Location','best');
title(ax1, sprintf('Total Cost (Z=%.3g)', Zl));
saveas(f1, fullfile(figDir,'costs.png'));

%% ---- 6) Map / Network (fallback line map) -------------------------------
f2 = figure('Color','w','Name','CCS Map'); ax2 = axes(f2); hold(ax2,'on');
pal = [0.30 0.30 0.30; 0.20 0.55 0.95; 0.10 0.70 0.40];
nodeColors = pal(max(cls,1),:);

tDraw = max(1,T);
F = qCO2_total(:,:,tDraw);
if any(F(:)>0)
    F(F<max(F(:))*1e-4) = 0;
    LW = 0.5 + 3.5 * (F / max(F(:)));
else
    LW = zeros(N,N)+0.5;
end

% Dashed if offshore (aco>1) or any maritime endpoint
isOff = (aco>1) | (mar_n(:)'>0.5) | (mar_n(:)>0.5);

for i=1:N
    for j=1:N
        if i==j || F(i,j)<=0, continue; end
        ls = '-'; if isOff(i,j), ls='--'; end
        hL = plot(ax2, [lon(i) lon(j)], [lat(i) lat(j)], ...
          'LineStyle', ls, 'LineWidth', LW(i,j), 'Color', [0 0 0]);  % RGB only
        % (Lines don't have per-object alpha in most MATLAB releases; leave as-is.)
    end
end
hSc = scatter(ax2, lon, lat, 28, nodeColors, 'filled', 'MarkerEdgeColor',[0 0 0]);  % RGB only
% Optional: edge transparency, if your MATLAB supports it (R2018a+)
if isprop(hSc,'MarkerEdgeAlpha'), set(hSc,'MarkerEdgeAlpha',0.4); end

xlabel(ax2,'Longitude'); ylabel(ax2,'Latitude');
title(ax2, sprintf('Network Flows (t=%d). Solid=Onshore, Dashed=Offshore', tDraw));
legend(ax2, {'Pipes','Nodes'}, 'Location','bestoutside');
xlim(ax2,[min(lon)-2,max(lon)+2]); ylim(ax2,[min(lat)-2,max(lat)+2]);
saveas(f2, fullfile(figDir,'map.png'));

%% ---- 7) Sankey / layered digraph ---------------------------------------
f3 = figure('Color','w','Name','Sankey'); ax3 = axes(f3); hold(ax3,'on');
[sIdx,tIdx,w] = local_arc_list(F);
if isempty(w)
    text(0.5,0.5,'No transport in selected year','HorizontalAlignment','center');
else
    if exist('SSankey','class') == 8
        S = SSankey(sIdx, tIdx, w); S.draw(ax3);
        title(ax3, sprintf('Sankey (t=%d)', tDraw));
    else
        G = digraph(sIdx, tIdx, w);
        p = plot(ax3, G, 'Layout','layered', 'ArrowSize', 10);
        p.LineWidth = 0.5 + 4*(G.Edges.Weight / max(G.Edges.Weight));
        title(ax3, sprintf('Flow graph (layered) — t=%d', tDraw));
    end
end
saveas(f3, fullfile(figDir,'sankey.png'));

fprintf('Figures saved to %s\n', figDir);

%% ============================ Local helpers ==============================
function sym = fetchSym(m,name)
% Try common accessors; stop if unavailable
if ismethod(m,'get'),       try sym = m.get(name);       return; end, end
if ismethod(m,'getSymbol'), try sym = m.getSymbol(name); return; end, end
if ismethod(m,'symbol'),    try sym = m.symbol(name);    return; end, end
if ismethod(m,'sym'),       try sym = m.sym(name);       return; end, end
error('Symbol "%s" not found via Transfer API.', name);
end

function T = trySymTable(m,name,cols)
try
    T = normalizeSymbolRecords(fetchSym(m,name), cols);
catch
    T = table();
end
end

function T = normalizeSymbolRecords(sym, wantCols)
% Convert a Transfer symbol's .records to a table with columns wantCols
if ~isprop(sym,'records'), T = table(); return; end
R = sym.records;
if istable(R)
    T = R;
elseif iscell(R)
    T = cell2table(R);
else
    T = table();
end
% Ensure column names exist (best-guess mapping)
if isempty(T), T = cell2table(cell(0,numel(wantCols))); end
n = width(T);
names = wantCols;
for k=1:min(n,numel(names))
    T.Properties.VariableNames{k} = names{k};
end
if n < numel(names)
    % pad missing columns
    for k=n+1:numel(names)
        T.(names{k}) = cell(height(T),1);
    end
end
% Ensure last column is named 'value'
T.Properties.VariableNames{end} = 'value';
end

function T = readCsvDump(gdxdumpCmd, gdxFileQuoted, sym, wantCols)
% Dumps a single symbol to CSV using gdxdump and reads it as a table.
tmp = [tempname '.csv'];
cmd = sprintf('%s %s symb=%s format=csv > "%s"', gdxdumpCmd, gdxFileQuoted, sym, tmp);
[status,outMsg] = system(cmd);
assert(status==0, 'gdxdump failed for "%s". Ensure GAMS on PATH. Message: %s', sym, outMsg);
opts = detectImportOptions(tmp,'NumHeaderLines',0);
T = readtable(tmp, opts);
delete(tmp);
% Normalize column names to wantCols length, last is 'value'
n = width(T);
for k=1:min(n,numel(wantCols))
    T.Properties.VariableNames{k} = wantCols{k};
end
if n < numel(wantCols)
    for k=n+1:numel(wantCols)
        T.(wantCols{k}) = cell(height(T),1);
    end
end
T.Properties.VariableNames{end} = 'value';
end

function exe = autodetect_gdxdump()
% Try common install locations, then PATH
candidates = {
    'C:\GAMS\gdxdump.exe'
    'C:\GAMS\48\gdxdump.exe'
    'C:\GAMS\49\gdxdump.exe'
    'C:\GAMS\50\gdxdump.exe'
};
exe = '';
for i=1:numel(candidates)
    if exist(candidates{i}, 'file'), exe = ['"' candidates{i} '"']; break; end
end
if isempty(exe)
    [status,out] = system('where gdxdump');
    if status==0
        lines = regexp(strtrim(out), '\r?\n', 'split');
        if ~isempty(lines) && exist(lines{1},'file'), exe = ['"' lines{1} '"']; end
    end
end
if isempty(exe)
    error(['No pude localizar gdxdump.exe.\n' ...
           'Solución rápida: establece su ruta manualmente o añade GAMS al PATH.']);
end
end

function s = toStrCol(x)
% Ensure column as cellstr
if iscell(x), s = cellfun(@char,x,'uni',false);
else, s = cellstr(string(x));
end
end

function M = tableToDense(R, domnames, domlabels, targetSz, allowEmpty)
% Build dense array from table R with (unknown order) domain columns + 'value'.
% domnames = canonical names (e.g., {'q','n','nn','t'})
if nargin<5, allowEmpty = false; end
if isempty(R)
    if allowEmpty, M = zeros(targetSz); return;
    else, error('Required symbol missing to build dense matrix.');
    end
end
VN = R.Properties.VariableNames;
% pick domain columns present (ignore 'value')
cand = intersect(domnames, VN, 'stable');
if numel(cand) ~= numel(domnames)
    % fallback: take first N non-value columns
    tmp = VN(~strcmpi(VN,'value'));
    cand = tmp(1:numel(domnames));
end
% label→index maps
maps = cell(1,numel(domnames));
for k=1:numel(domnames)
    labs = domlabels{k};
    key  = cellfun(@char, labs(:), 'uni', false);
    maps{k} = containers.Map(key, num2cell(1:numel(labs)));
end
% Match each cand column to canonical axis by label membership
col2axis = zeros(1,numel(domnames));
for k=1:numel(domnames)
    colk = cand{k};
    vals = R.(colk);
    if ~iscell(vals), vals = cellstr(string(vals)); end
    % pick first non-empty
    j = find(~cellfun(@isempty,vals), 1);
    if ~isempty(j)
        sample = vals{j};
        for a=1:numel(domnames)
            if isKey(maps{a}, sample), col2axis(k) = a; break; end
        end
    end
end
% Fill any unresolved by first-come assignment
missing = setdiff(1:numel(domnames), col2axis(col2axis>0));
col2axis(col2axis==0) = missing;

% Build dense
M = zeros(targetSz);
vname = 'value'; if ~ismember('value',VN), vname = VN{end}; end
H = height(R);
for i=1:H
    idx = cell(1,numel(domnames));
    skip = false;
    for k=1:numel(domnames)
        colk = cand{k};
        label = R.(colk){i};
        if isempty(label)
            skip = true; break;
        end
        if ~ischar(label), label = char(string(label)); end
        ax = col2axis(k);
            if ~isKey(maps{ax}, label)
            %   unseen label...
                skip = true; break;
            end
        idx{ax} = maps{ax}(label);
    end
    if skip, continue; end
    M(idx{:}) = R.(vname)(i);
end
% Final guard
szM = size(M); szM(end+1:numel(targetSz)) = 1;
if ~isequal(szM, targetSz)
    error('Dims mismatch after build: got %s, expected %s', mat2str(szM), mat2str(targetSz));
end
end

function v = vecFromTable(R, domnames, domlabels, targetSz)
v = tableToDense(R, domnames, domlabels, targetSz, true); v = v(:);
end

function x = scalarFromTable(R)
if isempty(R), x = NaN; return; end
vn = R.Properties.VariableNames;
if ismember('value', vn), x = R.value(1);
else, x = R{1,end};
end
end

function [sIdx,tIdx,w] = local_arc_list(F)
N = size(F,1); [sIdx,tIdx,w] = deal([]);
if ~any(F(:)>0), return; end
M = F>0 & ~eye(N);
[si,tj] = find(M);
w = F(sub2ind([N N],si,tj)); sIdx = si; tIdx = tj;
end
