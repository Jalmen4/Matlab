function [inlandRegionsMatrix, coastRegionsMatrix, maritimeRegionsMatrix, europeRegions, EU, EUbuf] = GettingLocationRegions( ...
    regionsCoordinates, variationLatitude, variationLongitude, latLim, lonLim, outputFormat)

% outputFormat: 'logical' (default) or 'numeric' for inland/coastal/maritime outputs.
% europeRegions is ALWAYS returned as logical.

if nargin < 6 || isempty(outputFormat)
    outputFormat = 'logical';
end

nRegionsDivisions = 20;

%% SUPPRESS WARNINGS
warning('off','MATLAB:polyshape:repairedBySimplify');
warning('off','map:shapefile:shxFileNotFound');
warning('off','map:shapefile:dbfFileNotFound');

%% SHAPEFILE (Natural Earth Admin 0)
shapefileFolder = 'C:\GIS\natural_earth\';
resolution = '50m';
shapefilePath = fullfile(shapefileFolder, ['ne_' resolution '_admin_0_countries.shp']);

nRowsRegionsCoordinates = size(regionsCoordinates,1);

%% STEP 1: Build subdivisions of regions (for classification robustness only)
regionsDivision = struct('Coordinates', cell(nRowsRegionsCoordinates,1));
for i = 1:nRowsRegionsCoordinates
    latVals = linspace(regionsCoordinates(i,1)-variationLatitude/2, ...
                       regionsCoordinates(i,1)+variationLatitude/2, nRegionsDivisions);
    lonVals = linspace(regionsCoordinates(i,2)-variationLongitude/2, ...
                       regionsCoordinates(i,2)+variationLongitude/2, nRegionsDivisions);

    [LAT, LON] = ndgrid(latVals, lonVals);
    coords = [LAT(:), LON(:)];
    insideBox = coords(:,1) >= latLim(1) & coords(:,1) <= latLim(2) & ...
                coords(:,2) >= lonLim(1) & coords(:,2) <= lonLim(2);
    regionsDivision(i).Coordinates = coords(insideBox,:);
end

%% STEP 2: Build Europe land poly (union of selected countries)
countries    = shaperead(shapefilePath, 'UseGeoCoords', true);
countryNames = regexprep({countries.ADMIN}, '\s+$', '');

europeCountryList = {
    'United Kingdom','Ireland','France','Spain','Portugal','Belgium','Netherlands','Germany',...
    'Denmark','Norway','Sweden','Finland','Estonia','Latvia','Lithuania','Poland','Czech','Slovak',...
    'Hungary','Austria','Switzerland','Italy','Slovenia','Croatia','Bosnia','Serbia',...
    'Montenegro','Kosovo','Albania','Macedonia','Greece','Bulgaria','Romania','Moldova',...
    'Ukraine','Belarus','Russia','Luxembourg','Andorra','Monaco','San Marino','Vatican','Liechtenstein','Turkey'
};

isEurope = false(size(countryNames));
for ii = 1:numel(europeCountryList)
    isEurope = isEurope | contains(countryNames, europeCountryList{ii}, 'IgnoreCase', true);
end
europeCountries = countries(isEurope);

% Polyshapes (clip Russia at 60E)
euPolys = polyshape.empty;
for k = 1:numel(europeCountries)
    x = europeCountries(k).Lon(:);
    y = europeCountries(k).Lat(:);
    nanBreaks = isnan(x) | isnan(y);
    startIdx = [1; find(nanBreaks)+1];
    endIdx   = [find(nanBreaks)-1; numel(x)];

    for part = 1:numel(startIdx)
        xi = x(startIdx(part):endIdx(part));
        yi = y(startIdx(part):endIdx(part));
        if numel(xi) < 3, continue; end

        if contains(countryNames{k},'Russia','IgnoreCase',true)
            keep = xi <= 60;
            xi = xi(keep); yi = yi(keep);
            if numel(xi) < 3, continue; end
            if xi(1)~=xi(end) || yi(1)~=yi(end)
                xi(end+1)=xi(1); yi(end+1)=yi(1);
            end
        end

        try
            p = polyshape(xi, yi, 'Simplify', false);
            if area(p) > 0
                euPolys(end+1) = p; %#ok<AGROW>
            end
        catch
        end
    end
end
assert(~isempty(euPolys), 'No Europe polygons built.');

EU   = simplify(polybuffer(union(euPolys),0), 'KeepCollinearPoints', true);

% Coastal buffer width ~ half the region radius
regionRadius = 0.5 * sqrt(variationLatitude^2 + variationLongitude^2);
coastalWidth = max(1e-6, 0.5 * regionRadius);
EUbuf = polybuffer(EU, coastalWidth);

%% STEP 3: Africa exclusion polygon
africaLon = [-5.60, -22.5, -22.50, 55, 55, 33, 22, 11.5, 11.5, 4.2, 1.17, 0.19, -0.15, -1.3];
africaLat = [35.95, 35.95, -37, -37, 15, 33, 34, 36, 38, 37.15, 36.75, 36.47, 36.27, 35.95];
[africaLon, ia] = unique(africaLon,'stable'); africaLat = africaLat(ia);
africaPoly = simplify(polyshape(africaLon, africaLat), 'KeepCollinearPoints', true);

%% STEP 4: Classification (using subdivision coverage + buffer)
Npts = nRegionsDivisions^2;
inlandFracThreshold   = 1 - (1 / Npts);     % almost all points on land
maritimeFracThreshold = min(4 / Npts, 0.05);% few land points allowed

% Build as logical internally; cast later if needed
inlandRegionsMatrix   = false(nRowsRegionsCoordinates,1);
coastRegionsMatrix    = false(nRowsRegionsCoordinates,1);
maritimeRegionsMatrix = false(nRowsRegionsCoordinates,1);
europeRegions         = false(nRowsRegionsCoordinates,1);

regionCenters = zeros(nRowsRegionsCoordinates,2);
regionClass   = strings(nRowsRegionsCoordinates,1);

for i = 1:nRowsRegionsCoordinates
    coords = regionsDivision(i).Coordinates;
    if isempty(coords), continue; end

    cLat = mean(coords(:,1));
    cLon = mean(coords(:,2));
    regionCenters(i,:) = [cLat, cLon];

    % Exclude Africa
    if isinterior(africaPoly, cLon, cLat)
        continue
    end
    europeRegions(i) = true;

    % Land overlap fraction
    inEU = isinterior(EU, coords(:,2), coords(:,1));
    fracInside = sum(inEU) / size(coords,1);

    if fracInside >= inlandFracThreshold
        inlandRegionsMatrix(i) = true;   regionClass(i) = "inland";
    elseif fracInside >= maritimeFracThreshold || isinterior(EUbuf, cLon, cLat)
        coastRegionsMatrix(i) = true;    regionClass(i) = "coastal";
    else
        maritimeRegionsMatrix(i) = true; regionClass(i) = "maritime";
    end
end

% Optional cast of class masks (europeRegions stays logical)
if strcmpi(outputFormat,'numeric')
    inlandRegionsMatrix   = double(inlandRegionsMatrix);
    coastRegionsMatrix    = double(coastRegionsMatrix);
    maritimeRegionsMatrix = double(maritimeRegionsMatrix);
end

%% STEP 5: Figures

% ===== FIGURE 1: solid coloured region blocks (one PATCH per region) =====
figure('Name','Europe Region Classification');
plot(EU, 'FaceColor',[0.9 0.9 0.9], 'FaceAlpha',0.5, 'EdgeColor','k'); hold on;

for i = 1:nRowsRegionsCoordinates
    if ~europeRegions(i), continue; end

    lat0 = regionsCoordinates(i,1);
    lon0 = regionsCoordinates(i,2);
    latMin = lat0 - variationLatitude/2;
    latMax = lat0 + variationLatitude/2;
    lonMin = lon0 - variationLongitude/2;
    lonMax = lon0 + variationLongitude/2;

    xi = [lonMin, lonMax, lonMax, lonMin];
    yi = [latMin, latMin, latMax, latMax];

    if inlandRegionsMatrix(i)
        patch(xi, yi, 'g', 'FaceAlpha',0.2, 'EdgeColor','k', 'LineWidth',0.25);
    elseif coastRegionsMatrix(i)
        patch(xi, yi, 'y', 'FaceAlpha',0.2, 'EdgeColor','k', 'LineWidth',0.25);
    elseif maritimeRegionsMatrix(i)
        patch(xi, yi, 'b', 'FaceAlpha',0.2, 'EdgeColor','k', 'LineWidth',0.25);
    end
end

% Legend (3 opaque squares)
hIn  = plot(nan, nan, 's','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',12);
hCo  = plot(nan, nan, 's','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',12);
hMar = plot(nan, nan, 's','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',12);

plot(EUbuf, 'FaceColor','none', 'EdgeColor',[0.3 0.3 0.3], 'LineStyle','--');
title('European Regions Classification');
xlabel('Longitude'); ylabel('Latitude');
axis([lonLim, latLim]); box on;
legend([hIn,hCo,hMar], {'Inland','Coastal','Maritime'}, 'Location','northwest','FontSize',11);
hold off;

% ===== FIGURE 2: diagnostic centers with BIG markers + sequential IDs =====
figure('Name','Diagnostic: Region Centers Classification');
plot(EU, 'FaceColor',[0.9 0.9 0.9], 'FaceAlpha',0.5, 'EdgeColor','k'); hold on;

markerSize = 300; % bigger so numbers fit
eIdx = find(europeRegions);            % indices of Europe-only regions
IDs  = 1:numel(eIdx);                  % sequential labels 1..Neurope

for kk = 1:numel(eIdx)
    i = eIdx(kk);
    if regionClass(i) == "inland"
        scatter(regionCenters(i,2), regionCenters(i,1), markerSize, 'o', ...
            'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
    elseif regionClass(i) == "coastal"
        scatter(regionCenters(i,2), regionCenters(i,1), markerSize, 'o', ...
            'MarkerFaceColor','y','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
    elseif regionClass(i) == "maritime"
        scatter(regionCenters(i,2), regionCenters(i,1), markerSize, 'o', ...
            'MarkerFaceColor','b','MarkerEdgeColor','k','MarkerFaceAlpha',0.2); % transparent blue
    end
    % Black label, centered
    text(regionCenters(i,2), regionCenters(i,1), sprintf('%d', IDs(kk)), ...
        'Color','k','FontSize',10,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

% Legend (exactly 3 items, using dummy handles)
hIn2  = scatter(nan, nan, markerSize, 'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
hCo2  = scatter(nan, nan, markerSize, 'o','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
hMar2 = scatter(nan, nan, markerSize, 'o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerFaceAlpha',0.2);

title('Diagnostic: Region Centers');
xlabel('Longitude'); ylabel('Latitude');
axis([lonLim, latLim]); grid on;
lgd = legend([hIn2,hCo2,hMar2], {'Inland','Coastal','Maritime'}, 'Location','northwest','FontSize',11);
lgd.ItemTokenSize = [20,20];
hold off;

end
