function createGDX(regionsEmission, regionsCoordinates, sequestrationCO2Capacity, ...
                   connectedRegionsMatrix, inlandRegionsMatrix, ...
                   coastRegionsMatrix, maritimeRegionsMatrix, ...
                   intraCellDist, regionsDistance, depthAquifers)

% Builds Model_Data.gdx for GAMS.
% UNITS (align with your GAMS model):
%   E, S, flows : [kt/year]
%   d           : [km]
%   UCC         : [€/t]
%   UTC         : [€/t·km]
%   USC         : [€/t]
%   aco, USC_offshore, cf, ce : [-] multipliers / shares

addpath("C:\GAMS\48\api\matlab")
import gams.transfer.*

%% --- 1) Sizes ---
nEuropeRegions = size(regionsEmission,1);

%% --- 2) Region names g1..gN ---
nRegions = arrayfun(@(i) sprintf('g%d',i), 1:nEuropeRegions, 'uni', false);

%% --- 3) Container & sets ---
m = Container();

t = Set(m, 't', 'records', {'t'}, 'description', 'years');
n = Set(m, 'n', 'records', nRegions, 'description', 'regions');
p = Set(m, 'p', 'records', {'Latitude','Longitude'}, 'description', 'coordinates');

% Cost / tech sets (must match GAMS exactly)
cNames = {'coalPostCombustion','gasPostCombustion','coalOxyCombustion','preCombustion'};
qNames = {'q1','q2','q3','q4'};
c = Set(m, 'c', 'records', cNames, 'description','capture technologies');
q = Set(m, 'q', 'records', qNames, 'description','flowrate discretisation');

%% --- 4) Parameters (core data) ---
cn        = Parameter(m,'cn',        {n,n}, 'description','Regions that are connected (0/1)');
coastal_n = Parameter(m,'coastal_n',  n,    'description','Coastal regions (0/1)');
mar_n     = Parameter(m,'mar_n',      n,    'description','Maritime regions (0/1)');
inland_n  = Parameter(m,'inland_n',   n,    'description','Inland regions (0/1)');

E     = Parameter(m,'E',     {n,t},  'description','Emissions [kt/y]');
coord = Parameter(m,'coord', {n,p},  'description','Coordinates [deg]');
S     = Parameter(m,'S',      n,     'description','Storage capacity [kt/y]');
depth = Parameter(m,'depth',  n,     'description','Storage depth [m]');
n_size= Parameter(m,'n_size', n,     'description','Region size [km or arbitrary]');
d     = Parameter(m,'d',     {n,n},  'description','Distance between regions [km]');

%% --- 5) Parameters (cost helpers & cost inputs) ---
aco          = Parameter(m,'aco',          {n,n}, 'description','Offshore pipeline factor [-]');
USC_offshore = Parameter(m,'USC_offshore',  n,    'description','Offshore sequestration multiplier [-]');

UCC = Parameter(m,'UCC', c, 'description','Unitary Capture Cost [€/t]');
UTC = Parameter(m,'UTC', q, 'description','Unit pipeline transport cost per distance [€/t·km]');
USC = Parameter(m,'USC', [], 'description','Unitary Sequestration Cost [€/t]');

%% --- 6) NEW: Capture feasibility & efficiency to export ---
cf  = Parameter(m,'cf', {n,c}, 'description','Tech feasibility in region n (share) [-]');
ce  = Parameter(m,'ce', c,     'description','Capture efficiency of tech c [-]');

%% --- 7) Ensure shape (column vectors) for masks ---
coastRegionsMatrix    = double(coastRegionsMatrix(:));
maritimeRegionsMatrix = double(maritimeRegionsMatrix(:));
inlandRegionsMatrix   = double(inlandRegionsMatrix(:));

%% --- 8) Build USC_offshore and aco from masks ---
offshoreSequestrationMultiplier = 2.5;  % USC_offshore = 2.5 if maritime, else 1
offshorePipelineFactor          = 1.5;  % aco(i,j) = 1.5 if any endpoint maritime; else 1

isMar = maritimeRegionsMatrix == 1;

USC_offshore_vec = ones(nEuropeRegions,1);
USC_offshore_vec(isMar) = offshoreSequestrationMultiplier;

acoMat = ones(nEuropeRegions, nEuropeRegions);
for i = 1:nEuropeRegions
    for j = 1:nEuropeRegions
        if i ~= j && connectedRegionsMatrix(i,j) == 1
            if isMar(i) || isMar(j)
                acoMat(i,j) = offshorePipelineFactor;   % maritime endpoint → offshore markup
            else
                acoMat(i,j) = 1;                        % coastal–coastal & inland combos → onshore
            end
        end
    end
end

%% --- 9) Cost values ---
% UCC(c) [€/t] (order follows cNames)
uccVals = [33; 54; 36; 25];

% UTC(q) [€/t·km] (onshore pipeline by level; order follows qNames)
utcVals = [0.054; 0.028; 0.016; 0.010];

% USC [€/t]
uscScalar = 7.2;

%% --- 10) NEW: ce(c) and cf(n,c) values ---
% ce(c) [-] (order follows cNames)
ceVals = [0.87; 0.88; 0.92; 0.86];

% cf(n,c) [-] set uniformly to 0.5 for every (n,c)
Ni = nEuropeRegions;
Nc = numel(cNames);
[nIdx, cIdx] = ndgrid(1:Ni, 1:Nc);
cf_n_keys = categorical(nRegions(nIdx(:)));
cf_c_keys = categorical(cNames(cIdx(:)));
cf_vals    = 0.5 * ones(Ni * Nc, 1);

%% --- 11) Assign values to all parameters ---
cn.setRecords(connectedRegionsMatrix);
coastal_n.setRecords(coastRegionsMatrix);
mar_n.setRecords(maritimeRegionsMatrix);
inland_n.setRecords(inlandRegionsMatrix);

E.setRecords(regionsEmission);
coord.setRecords(regionsCoordinates);
S.setRecords(sequestrationCO2Capacity);
depth.setRecords(depthAquifers);
n_size.setRecords(intraCellDist);
d.setRecords(regionsDistance);

aco.setRecords(acoMat);
USC_offshore.setRecords(USC_offshore_vec);

% 1D parameters: struct with domain (categorical column) + value
UCC.setRecords( struct('c', categorical(cNames(:)), 'value', uccVals(:)) );
UTC.setRecords( struct('q', categorical(qNames(:)), 'value', utcVals(:)) );
USC.setRecords( uscScalar );

ce.setRecords(  struct('c', categorical(cNames(:)), 'value', ceVals(:)) );

% 2D parameter: struct with both domain columns + value (all rows)
cf.setRecords( struct('n', cf_n_keys, 'c', cf_c_keys, 'value', cf_vals) );

%% --- 12) Write GDX ---
m.write('Model_Data.gdx');

end
