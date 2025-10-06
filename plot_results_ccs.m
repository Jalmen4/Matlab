%% plot_results_ccs_script.m
% Self-contained script to visualize CCS results from Results.gdx
% (1) Cost breakdown
% (2) Europe lon–lat cells + transport network (on/offshore)
% (3) Sankey of capture → transport → storage (uses SSankey if available, digraph fallback)

clc;

%% --- Paths (adjust if needed)
addpath('C:\GAMS\48\api\matlab');     % GAMS Transfer API
% addpath('C:\path\to\SSankey');      % If you use SSankey class

%% --- Inputs
resultsGdx = 'Results.gdx';            % ensure this exists in current folder
timeIdx    = 1;                        % choose year index to plot

% If you prefer loading geo from file, uncomment:
% load('geo_inputs.mat', ...
%   'regionsCoordinates','variationLatitude','variationLongitude','latLim','lonLim', ...
%   'inlandRegionsMatrix','coastRegionsMatrix','maritimeRegionsMatrix','europeRegions', ...
%   'EU','EUbuf');

% Geo struct from variables already in workspace
geo = struct();
geo.regionsCoordinates    = regionsCoordinates;
geo.variationLatitude     = variationLatitude;
geo.variationLongitude    = variationLongitude;
geo.latLim                = latLim;
geo.lonLim                = lonLim;
geo.inlandRegionsMatrix   = inlandRegionsMatrix(:);
geo.coastRegionsMatrix    = coastRegionsMatrix(:);
geo.maritimeRegionsMatrix = maritimeRegionsMatrix(:);
geo.europeRegions         = europeRegions(:);
if exist('EU','var'),    geo.EU    = EU;    end
if exist('EUbuf','var'), geo.EUbuf = EUbuf; end

%% --- Read GDX via GAMS Transfer
m = gams.transfer.Container();
m.read(resultsGdx, 'format', 'dense_matrix');
S = m.data;

%% --- Labels & core data
Nlab_raw = string(S.n.getUELs(1));
Clab     = string(S.c.getUELs(1));
Plab     = string(S.p.getUELs(1));
Tlab     = string(S.t.getUELs(1));

coord = S.coord.records.value;                 % (#N x #p) [Latitude, Longitude]
d     = S.d.records.value;                     % (#N x #N)
aco   = S.aco.records.value;                   % (#N x #N)
USC_offshore = S.USC_offshore.records.value(:);
UCC   = S.UCC.records.value(:);
UTC   = S.UTC.records.value(:);
USC   = S.USC.records.value;

% Results (levels exported as parameters)
TCC = S.TCC_l.records.value(:);
TTC = S.TTC_l.records.value(:);
TSC = S.TSC_l.records.value(:);
TC  = S.TC_l.records.value(:);
Z   = S.Z_l.records.value;

cCO2 = S.cCO2_l.records.value;                 % (#c x #N x #t) [kt/y] (axis order may vary)
qCO2 = S.qCO2_l.records.value;                 % unknown axis order
Seq  = S.Seq_l.records.value;                  % (#N x #t) [kt/y]

% Coordinates from GDX (master node order)
latIdx = find(strcmpi(Plab,'Latitude'),1);
lonIdx = find(strcmpi(Plab,'Longitude'),1);
assert(~isempty(latIdx) && ~isempty(lonIdx), 'Set p must contain Latitude & Longitude.');
latG = coord(:,latIdx);
lonG = coord(:,lonIdx);
N    = size(coord,1);

% Robust region labels in case UELs < or > N
if numel(Nlab_raw) < N
    pad = "g" + string((numel(Nlab_raw)+1):N);
    Nlab = [Nlab_raw(:); pad(:)];
elseif numel(Nlab_raw) > N
    Nlab = Nlab_raw(1:N);
else
    Nlab = Nlab_raw;
end

%% --- GEO reconciliation (map masks to GDX node order)
rc_geo   = geo.regionsCoordinates;   % [Ng x 2] (lat, lon)
vLat     = geo.variationLatitude;    % deg
vLon     = geo.variationLongitude;   % deg
latLim   = geo.latLim;
lonLim   = geo.lonLim;
EU       = []; if isfield(geo,'EU'),    EU    = geo.EU;    end
EUbuf    = []; if isfield(geo,'EUbuf'), EUbuf = geo.EUbuf; end

inland_g   = logical(geo.inlandRegionsMatrix(:));
coastal_g  = logical(geo.coastRegionsMatrix(:));
maritime_g = logical(geo.maritimeRegionsMatrix(:));
europe_g   = logical(geo.europeRegions(:));

Ng = size(rc_geo,1);
fprintf('[plot_results_script] N (GDX)=%d, Ng (geo)=%d\n', N, Ng);

rcGDX = [latG, lonG];
rcGeo = [rc_geo(:,1), rc_geo(:,2)];

idxGeo = zeros(N,1);
useIdentity = false;
if Ng == N
    dLat0 = abs(rcGeo(:,1) - rcGDX(:,1));
    dLon0 = abs(rcGeo(:,2) - rcGDX(:,2));
    maxDL = max(dLat0(:)); maxDN = max(dLon0(:));
    fprintf('[plot_results_script] max |Δlat|=%.6f°, max |Δlon|=%.6f°\n', maxDL, maxDN);
    tol = 1e-3; % ~110 m
    if maxDL <= tol && maxDN <= tol
        useIdentity = true;
        idxGeo = (1:N).';
        fprintf('[plot_results_script] Using 1:1 alignment.\n');
    end
end
if ~useIdentity
    fprintf('[plot_results_script] Using nearest-neighbor mapping GDX→GEO.\n');
    if exist('knnsearch','file') == 2
        idxGeo = knnsearch(rcGeo, rcGDX, 'K', 1);
    else
        for k = 1:N
            d2 = (rcGeo(:,1)-rcGDX(k,1)).^2 + (rcGeo(:,2)-rcGDX(k,2)).^2;
            [~, idxGeo(k)] = min(d2);
        end
    end
    dLat = abs(rcGeo(idxGeo,1) - rcGDX(:,1));
    dLon = abs(rcGeo(idxGeo,2) - rcGDX(:,2));
    if max(dLat) > max(vLat, 0.6) || max(dLon) > max(vLon, 0.6)
        warning('Geo–GDX mismatch looks large (max |Δlat|=%.3f°, |Δlon|=%.3f°).', max(dLat), max(dLon));
    end
end

% Masks in GDX order
inland   = inland_g(idxGeo);
coastal  = coastal_g(idxGeo);
maritime = maritime_g(idxGeo);
inEU     = europe_g(idxGeo);
if ~any(inEU)
    warning('inEU mask empty; drawing all nodes.');
    inEU = true(N,1);
end

% Draw coordinates
lat = latG; lon = lonG;

fprintf('[plot_results_script] inEU true count: %d of %d\n', sum(inEU), numel(inEU));

%% --- Output folder
figDir = fullfile(pwd,'figures'); if ~exist(figDir,'dir'), mkdir(figDir); end

%% ========== Figure 1: Cost breakdown ==========
f1 = figure('Color','w','Position',[100 100 920 440]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (a) Stacked bars
nexttile;
X = (1:numel(TCC))';  B = [TCC TTC TSC];
bar(X, B, 'stacked','LineWidth',0.6);
xlabel('Time period'); ylabel('Cost [€ / year]');
title('Total Cost Breakdown by Year');
legend({'Capture','Transport','Sequestration'},'Location','northoutside','Orientation','horizontal','Box','off');
box off; grid on; set(gca,'XTick',X,'XTickLabel',Tlab,'FontName','Times','FontSize',11);

% (b) Total cost line
nexttile;
plot(X, TC, 'LineWidth',1.4); grid on; box off;
xlabel('Time period'); ylabel('Total Cost [€ / year]');
title(sprintf('Total Cost per Year (Z = %.3g €)', Z));
set(gca,'XTick',X,'XTickLabel',Tlab,'FontName','Times','FontSize',11);

exportgraphics(f1, fullfile(figDir,'cost_breakdown.png'), 'Resolution',300);
exportgraphics(f1, fullfile(figDir,'cost_breakdown.pdf'));

%% ========== Figure 2: Map with cells + network ==========
% Robust aggregation of transport across q for selected time
Q = numel(UTC);
T = numel(Tlab);

qraw = qCO2;
totalQ = full(sum(qraw(:)));
fprintf('[plot_results_script] qCO2 total sum (all dims): %.6g kt/y\n', totalQ);

F = local_sum_over_q_select_time(qraw, N, T, Q, timeIdx);  % [N x N]

Active      = F > 1e-9;
OffshoreArc = (aco > 1.0001);
[iIdx, jIdx] = find(Active);
flows = F(Active);
isOff  = OffshoreArc(Active);
fprintf('[plot_results_script] Active arcs: %d | Offshore among active: %d\n', numel(flows), nnz(isOff));

% Edge widths
if ~isempty(flows)
    lw = 0.5 + 3.0 * (flows - min(flows)) / max(eps, (max(flows)-min(flows)));
else
    lw = [];
end

f2 = figure('Color','w','Position',[80 80 1100 820]); hold on;

% (a) EU polygons
if exist('EU','var') && ~isempty(EU)
    plot(EU,   'FaceColor',[0.92 0.92 0.92], 'FaceAlpha',0.7, ...
              'EdgeColor',[0.3 0.3 0.3], 'LineWidth',0.7);
end
if exist('EUbuf','var') && ~isempty(EUbuf)
    plot(EUbuf,'FaceColor','none', 'EdgeColor',[0.4 0.4 0.4], ...
               'LineStyle','--', 'LineWidth',0.6);
end

% (b) Cells: faint all, then stronger for inEU
colIn  = [0.10 0.65 0.10];
colCo  = [0.95 0.85 0.15];
colMar = [0.25 0.45 0.95];
edgeA_faint  = 0.25; edgeA_strong = 0.75;
lineW_faint  = 0.3;  lineW_strong = 0.6;

% all faint
for i = 1:N
    lat0 = lat(i); lon0 = lon(i);
    vlat = [lat0 - vLat/2, lat0 - vLat/2, lat0 + vLat/2, lat0 + vLat/2, lat0 - vLat/2];
    vlon = [lon0 - vLon/2, lon0 + vLon/2, lon0 + vLon/2, lon0 - vLon/2, lon0 - vLon/2];
    if inland(i),        cEdge = colIn;
    elseif coastal(i),   cEdge = colCo;
    elseif maritime(i),  cEdge = colMar;
    else,                cEdge = [0.5 0.5 0.5];
    end
    plot(vlon, vlat, '-', 'Color', [cEdge edgeA_faint], 'LineWidth', lineW_faint);
end

% inEU stronger
for i = find(inEU(:))'
    lat0 = lat(i); lon0 = lon(i);
    vlat = [lat0 - vLat/2, lat0 - vLat/2, lat0 + vLat/2, lat0 + vLat/2, lat0 - vLat/2];
    vlon = [lon0 - vLon/2, lon0 + vLon/2, lon0 + vLon/2, lon0 - vLon/2, lon0 - vLon/2];
    if inland(i),        cEdge = colIn;
    elseif coastal(i),   cEdge = colCo;
    elseif maritime(i),  cEdge = colMar;
    else,                cEdge = [0.5 0.5 0.5];
    end
    plot(vlon, vlat, '-', 'Color', [cEdge edgeA_strong], 'LineWidth', lineW_strong);
end

% (c) Pipes (onshore/offshore)
hOn  = plot(nan,nan,'-','Color',[0 0 0 0.28],'LineWidth',1.5); % legend proxy
hOff = plot(nan,nan,'--','Color',[0 0 0 0.28],'LineWidth',1.5);
for k = 1:numel(iIdx)
    i = iIdx(k); j = jIdx(k);
    if ~(inEU(i) && inEU(j)), continue; end
    style = '-'; if isOff(k), style = '--'; end
    plot([lon(i) lon(j)], [lat(i) lat(j)], style, ...
         'Color', [0 0 0 0.35], 'LineWidth', lw(max(k,1)));
end

% (d) Nodes & labels
cap_t  = squeeze(sum(cCO2(:,:,timeIdx), 1));    % [N x 1] kt/y (assumes c first)
if isrow(cap_t), cap_t = cap_t.'; end
inflow = sum(F, 1)';                            % [N x 1] kt/y
sizeBase = cap_t + inflow;
if all(sizeBase==0), sizeBase = ones(N,1); end
nodeSize = 40 + 60 * sqrt(sizeBase / max(1e-9, max(sizeBase)));
scatter(lon(inEU), lat(inEU), max(nodeSize(inEU),20), ...
        'filled', 'MarkerFaceAlpha',0.85, 'MarkerFaceColor',[0.2 0.2 0.2]);

idxEU = find(inEU);
for h = 1:numel(idxEU)
    ii = idxEU(h);
    text(lon(ii), lat(ii), " "+Nlab(ii), 'FontSize', 7, 'Color', [0.1 0.1 0.1]);
end

axis([lonLim latLim]); axis xy; grid on; box on;
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
title(sprintf('CO_2 Flow Network with Region Cells (time = %s)', Tlab(timeIdx)), ...
      'FontName','Times','FontSize',13);

legend([hOn hOff], {'Onshore','Offshore'}, ...
       'Location','southoutside','Orientation','horizontal','Box','off');

exportgraphics(f2, fullfile(figDir, sprintf('map_cells_t%d.png', timeIdx)), 'Resolution',300);
exportgraphics(f2, fullfile(figDir, sprintf('map_cells_t%d.pdf',  timeIdx)));

%% ========== Figure 3: Sankey (SSankey or digraph fallback) ==========
% Capture by tech (kt/y)
Ccap_c = squeeze(sum(cCO2(:,:,timeIdx), 2));   % [#c x 1]

% Transport split (on/offshore)
TotalFlow   = sum(F(:));
OffFlow     = sum(F(aco > 1.0001));
OnFlow      = TotalFlow - OffFlow;
shareOff = local_safeDiv(OffFlow, TotalFlow);
shareOn  = 1 - shareOff;

% Storage split (on/offshore)
Seq_t   = Seq(:, timeIdx);
OffStor = sum(Seq_t(USC_offshore > 1+1e-9));
OnStor  = sum(Seq_t(USC_offshore <= 1+1e-9));

% Nodes: techs + Trans(On) + Trans(Off) + Stor(On) + Stor(Off)
ClabRow    = reshape(Clab, 1, []);
nodeLabels = [ClabRow, "Transport (Onshore)", "Transport (Offshore)", "Storage (Onshore)", "Storage (Offshore)"];
Nc = numel(ClabRow);
TR_ON  = Nc+1;  TR_OFF = Nc+2;  ST_ON = Nc+3; ST_OFF = Nc+4;

% Adjacency
BN = Nc + 4;
Adj = zeros(BN, BN);
for i = 1:Nc
    cap_i = Ccap_c(i);
    if cap_i > 0
        Adj(i, TR_ON)  = Adj(i, TR_ON)  + cap_i*shareOn;
        Adj(i, TR_OFF) = Adj(i, TR_OFF) + cap_i*shareOff;
    end
end
if OnStor  > 0, Adj(TR_ON,  ST_ON)  = OnStor;  end
if OffStor > 0, Adj(TR_OFF, ST_OFF) = OffStor; end
Layer = [ones(1,Nc), 2, 2, 3, 3];

f3 = figure('Color','w','Position',[100 100 1000 520]);
ax3 = axes('Parent', f3); hold(ax3, 'on');
try
    % Try SSankey if available
    SK = SSankey(ax3, [], [], [], 'NodeList', cellstr(nodeLabels), 'AdjMat', Adj);
    SK.Layer              = Layer;
    SK.LabelLocation      = 'right';
    SK.ValueLabelLocation = 'none';
    SK.BlockScale         = 0.08;
    SK.Sep                = 0.06;
    SK.Align              = 'center';
    SK.RenderingMethod    = 'interp';
    SK.ColorList          = repmat([0.18 0.18 0.18], BN, 1);
    SK.draw();
    title(sprintf('Capture \\rightarrow Transport \\rightarrow Storage (time = %s)', Tlab(timeIdx)), ...
          'FontName','Times','FontSize',13);
catch ME
    warning('SSankey failed: %s. Falling back to layered digraph.', ME.message);
    [si, tj, val] = find(Adj);
    if isempty(val)
        text(ax3, 0.5, 0.5, 'No flows to display', 'Units','normalized', ...
            'HorizontalAlignment','center','FontSize',12);
    else
        G = digraph(si, tj, val, BN);
        p = plot(G, 'Layout','layered', 'Interpreter','none', 'Parent', ax3);
        % Robust linewidth scaling
        LW = 1.5;
        if istable(G.Edges)
            vn = string(G.Edges.Properties.VariableNames);
            wcol = find(strcmpi(vn,'Weight') | strcmpi(vn,'Weight_'), 1, 'first');
            if ~isempty(wcol)
                W = G.Edges.(vn(wcol));
                LW = 0.75 + 4*(W - min(W))/max(eps,(max(W)-min(W)));
            end
        end
        p.LineWidth = LW;
        p.NodeLabel = cellstr(nodeLabels(:));
    end
    title(sprintf('Capture → Transport → Storage (time = %s) [graph view]', Tlab(timeIdx)), ...
          'FontName','Times','FontSize',13);
end

exportgraphics(f3, fullfile(figDir, sprintf('sankey_t%d.png', timeIdx)), 'Resolution',300);
exportgraphics(f3, fullfile(figDir, sprintf('sankey_t%d.pdf',  timeIdx)));

fprintf('Saved figures in %s\n', figDir);

%% ================== Local helpers ==================
function r = local_safeDiv(a,b)
    if b <= 0, r = 0; else, r = a/b; end
end

function F = local_sum_over_q_select_time(qraw, N, T, Q, timeIdx)
    % Coerce qCO2 to [Q x N x N x T], then sum over Q and slice time
    q4 = local_coerce_qnnnt(qraw, N, T, Q);
    F_allT = squeeze(sum(q4, 1));    % [N x N x T] (or [N x N] if T=1)
    if ndims(F_allT) == 2
        F = F_allT;
    else
        F = F_allT(:,:,timeIdx);
    end
end

function q4 = local_coerce_qnnnt(qraw, N, T, Q)
    dims = size(qraw);
    if numel(dims) < 4, dims(end+1:4) = 1; end
    idxN  = find(dims == N);
    idxT  = find(dims == T, 1, 'first');
    idxQ  = find(dims == Q, 1, 'first');
    if isempty(idxQ)
        cand = setdiff(1:4, [idxN(:).' idxT]);
        if ~isempty(cand)
            szcand = dims(cand);
            mask = (szcand ~= 1) & (szcand ~= N) & (szcand ~= T);
            if any(mask)
                pos = find(mask, 1, 'first');
                idxQ = cand(pos);
            else
                idxQ = cand(1);
            end
        else
            idxQ = 1;
        end
    end
    if numel(idxN) < 2 || isempty(idxT) || isempty(idxQ)
        warning('qCO2 axes ambiguous. dims=%s (N=%d,T=%d,Q=%d). Best guess.', mat2str(dims), N, T, Q);
        need = [Q N N T];
        perm = zeros(1,4);
        left = 1:4;
        for k = 1:4
            [~,j] = min(abs(dims(left)-need(k)));
            perm(k) = left(j);
            left(j) = [];
        end
        q4 = permute(qraw, perm);
        return
    end
    if numel(idxN) > 2, idxN = idxN(1:2); end
    perm = [idxQ idxN(1) idxN(2) idxT];
    q4 = permute(qraw, perm);
end
