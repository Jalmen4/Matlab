%% SCRIPT USED TO BE ABLE TO CHANGE THE RESOLUTION OF THE REGIONS CONSIDERED AND CONSEQUENTLY MODIFY DE VALUES OF THE PARAMETERS USED IN THE MODEL (Regions CO2 emissions, Regions CO2capacity, surrounding regions, regions locations, and so on). 
% This is the main script. The other functions which are call from this
% main script are the next one:
% - ExtractingStorageData: Used to obtain the storage capacity of each
% region considered. The data are gotten from a official database.
% - utm2deg: Used to convert the coordinates obtained using "ExtractingStorageData" from UTM to geographical coordinates.
% - GettingLocationRegions: Used to find out the location of each region (Inland, coastal or maritime).
% - surroundingRegionsScript: Used to identify the regions that are surrounding another one.
% - regionsAndDataGDX: Used to create a gdx file will all the information.


clc
clear
close all

%% SETTING THE REGION CONSIDERED LIMITS AND THE RESOLUTION DESIRED

tic

nLatitudSlices = 11; % Division made in the latitude direction
nLongitudSlices = 11; % Division made in the longitude direction

% southLatitude = 35.95 ; %South Latitude Limit (degrees)
% northLatitude = 44 ;  %North Latitude Limit (degrees)
% 
% westLongitude = -9.75 ; % West Longitude Limit (degrees)
% eastLongitude = 4.5 ; % East Longitude Limit (degrees)

southLatitude = 35.95 ; %South Latitude Limit (degrees)
northLatitude = 70 ;  %North Latitude Limit (degrees)

westLongitude = -10; % West Longitude Limit (degrees)
eastLongitude = 30 ; % East Longitude Limit (degrees)

latLim = [southLatitude, northLatitude] ; % Latitude limit (South and North limits)
lonLim = [westLongitude, eastLongitude] ; % Longitude limit (West and East limits)


%% GETTING COORDINATES AND EMISSIONS DATA FROM EDGAR DATABASE 
% (https://edgar.jrc.ec.europa.eu/dataset_ghg70)

%folderNCFile = 'C:\Users\jam10\OneDrive - UNIVERSIDAD ALICANTE\Escritorio\OneDrive - UNIVERSIDAD ALICANTE\Doctorado\Tesis\MATLAB\Regions resolution increasing\ObtainingEmissionsLast50years\Emissions Data\' ; % Personal computer path
folderNCFile = 'C:\Users\josea\OneDrive - UNIVERSIDAD ALICANTE\Escritorio\OneDrive - UNIVERSIDAD ALICANTE\Doctorado\Tesis\MATLAB\Regions resolution increasing\ObtainingEmissionsLast50years\Emissions Data Flux\' ; % University computer path


ncFile = "v8.0_FT2022_GHG_CO2_2022_TOTALS_flx.nc" ; % Opened NC file
%ncFile = "v7.0_FT2021_CO2_excl_short-cycle_org_C_2021_TOTALS.0.1x0.1.nc" ; % Opened NC file

% Complete Path
fullFilename = strcat(folderNCFile,ncFile) ;

emissionsLongitudeNC = ncread(fullFilename,"lon")  ; % Degrees (0º - 360º), Size (3600 x 1)
emissionsLatitudeNC = ncread(fullFilename,"lat") ; % Degrees (-90º - 90º), Size (1800 x 1)
emissionsValueNC = ncread(fullFilename,"fluxes") ; % kg/(s *m2) being m2 the area of the grid considered (0.1º x 0.1º), Size (3600 x 1800) 

%% CHANGING LONGITUDE RANGE FROM (0º - +360º) TO (-180º - +180ºC) 

% The range of the Longitude data got from Edgar are from 0 to 360ºC, but as 
% usually longitude are expresed in a range from -180º to +180, its units
% are going to be changed for comfort

% Obtaining the rows and column numbers  of the Longitude and Latitude matrices
[emissionsLongitudeNCRows, emissionsLongitudeNCColumns] = size(emissionsLongitudeNC) ;
[emissionsLatitudeNCRows, emissionsLatitudeNCColumns] = size(emissionsLatitudeNC) ;
[emissionsValueNCRows, emissionsValueNCColumns] = size(emissionsValueNC) ;



for i=1:emissionsLongitudeNCRows
    if emissionsLongitudeNC(i,1)>=180
        emissionsLongitudeNC(i,1) = emissionsLongitudeNC(i,1) - 360 ;
    end
end



%% DIMINISHING DATA SIZE IN ORDER TO REDUCE THE RUNNING TIME 

% It is going to be done some modifications in the data obtained from EDGAR
% with the purpose of just taking into account the values belonging to the
% region considered by using its latitude and longitude limits.

% The new variables that are going to include the values of the considered region, 
% are going to be created in order to save time

emissionsLongitude = zeros(emissionsLongitudeNCRows,emissionsLongitudeNCColumns) ;
emissionsLatitude = zeros(emissionsLatitudeNCRows,emissionsLatitudeNCColumns) ;
emissionsValue = zeros(emissionsLongitudeNCRows,emissionsLatitudeNCRows) ;


% Loops are going to be used to build the new matrices. The function of
% these loops is to give a value of -1000 to the elements of the matrix that
% are outside the region (For the case of the Longitude & Latitude
% matrices), and a value of -1 in the case of the emission data.

for i=1:emissionsLongitudeNCRows
        if  emissionsLongitudeNC(i,:) < lonLim(1) || emissionsLongitudeNC(i,:) > lonLim(2)
            emissionsLongitude(i,:) = -1000 ;
            emissionsValue(i,:) = -1 ;

        else
            emissionsLongitude(i,:) = emissionsLongitudeNC(i,:) ;
            emissionsValue(i,:) = emissionsValueNC(i,:) ;
        end
end

for i=1:emissionsLatitudeNCRows
        if  emissionsLatitudeNC(i,:) < latLim(1) || emissionsLatitudeNC(i,:) > latLim(2)
            emissionsLatitude(i,:) = -1000 ;
            emissionsValue(:,i) = -1 ;

        else
            emissionsLatitude(i,:) = emissionsLatitudeNC(i,:) ;
            emissionsValue(:,i) = emissionsValue(:,i) ;
        end
end

% Once the values of the regions have been totally isolated, the next loop is used to get from the matrix just
% those values.

k = 1 ;

for i=1:emissionsValueNCRows 
   if all(emissionsValue(i,:) == -1) == false % Condition: If not all the values of the row i are equal to -1 then...
       regionPosition = find(emissionsValue(i,:) > -1) ; % Saving the column number of the values higher than -1 in a new vector called regionPosition
       emissionsValue2(k,:) = emissionsValue(i,regionPosition) ; % Creating the new matrix with these values
       k = k + 1 ;
   end

end

clear emissionsValue emissionsLongitudeNC emissionsLatitudeNC emissionsValueNC

emissionsValue = emissionsValue2 ;

% The same that has been done with the emissionValue matrix need to be
% carried out with the longitude and latitude matrices, but it is not
% necessary to employ loops.

LongitudePosition = find(emissionsLongitude(:,1) > -1000) ;
emissionsLongitude = emissionsLongitude(LongitudePosition,1) ;

LatitudePosition = find(emissionsLatitude(:,1) > -1000) ;
emissionsLatitude = emissionsLatitude(LatitudePosition,1) ;




[emissionsLongitudeRows, emissionsLongitudeColumns] = size(emissionsLongitude) ;
[emissionsLatitudeRows, emissionsLatitudeColumns] = size(emissionsLatitude) ;



%% GETTING A COORDINATES MATRIX 

% Obtaining a coordinates matrix based on the Longitude and Latitude
% matrices that includes all the possible combinations

it = 0 ;

% In order to reduce the resolution time, it is prefered to create the
% Coordinates matrix previously (It is possible to know its size
% beforehand because the size of its origin matrices it is known {Longitude &
% Latitude}).
gridCoordinates = zeros(emissionsLongitudeRows * emissionsLatitudeRows, emissionsLongitudeColumns + emissionsLatitudeColumns) ; 

for i=1:emissionsLatitudeRows
    for j=1:emissionsLongitudeRows
        it = it + 1 ;
        gridCoordinates(it,:) = [emissionsLatitude(i,1), emissionsLongitude(j,1)] ;
    end
end

[gridCoordinatesRows, gridCoordinatesColumns] = size(gridCoordinates) ;





%% GRID AREA (0.1º X 0.1º) EXPRESSED INTO SQUARE METERS

% To express the emission data in mass/time units, it is needed to multiply
% the values obtained from EDGAR database by the area of the grid considered.


% Knowing that the coordinates are referred to the center of each
% grid-cell, the area is going to be calculated as the area of a square.
% In order to do that it is needed to establish 3 apexes of each grid:
% ARefPoint(lower left corner)
% BRefPoint(upper left corner)
% CRefPoint(lower right corner)

refDifValue = 0.05 ; % Units (Degrees)
ARefPoint = gridCoordinates - refDifValue ;
BRefPoint = [ARefPoint(:,1) + 0.1, ARefPoint(:,2)] ;
CRefPoint = [ARefPoint(:,1), ARefPoint(:,2) + 0.1] ;


% Once the grid apexes are defined, it is possible to calculate the area
% based on the sides length (Distance between the reference points)

wgs84 = wgs84Ellipsoid("m") ; % Line to create a World Geodetic reference Ellipsoid  with length in meters.
distBA = distance(BRefPoint,ARefPoint,wgs84) ;
distCA = distance(CRefPoint,ARefPoint,wgs84) ;
gridArea = distBA.*distCA ; % m2

clear CRefPoint BRefPoint ARefPoint distBA distCA



%% CREATING A MATRIX (dataEDGAR) WHICH INCLUDES ALL THE CELLS AND ITS RESPECTIVE EMISSIONS

dataEdgar =  zeros(gridCoordinatesRows, gridCoordinatesColumns + 1) ;

% Knowing the size of the gridCoordinate and emissionsValue, it is possible
% to know how to build the matrixDataEdgar.

dataEdgarSlices = gridCoordinatesRows/emissionsLongitudeRows ;



y = 1;
x = 0 ;
for i = 1:dataEdgarSlices
    x= x + 1 ;
    a = emissionsLongitudeRows ;
    dataEdgar(y:y+a-1,:) = [gridCoordinates(y:y+a-1,:), emissionsValue(:,x)] ;

    y = y + a ;
   

end


%% OBTAINED DATA EDGAR MATRIX WITH EMISSIONS DATA EXPRESSED IN TON OF CO2/YEAR

SECOND2YEAR = 60*60*24*365 ; % Second to year conversion factor Units (s/year)
dataEdgar = [dataEdgar(:,1:2), dataEdgar(:,3).*gridArea.*SECOND2YEAR./1000] ; % Conversion of emission elements from kg/(s*m2) to ton/year. Units [Degrees,Degrees,tonCO2/year]

totalCO2 = 3.77521e+10; % Total CO2 value that is showns in the EDGAR file. Units (ton/year)
totalCO2Calculated = sum(dataEdgar(:,3)) ; % Total CO2 value obtained using this program. Units (ton/year)


%% BUILDING THE REGIONS COORDINATES MATRIX

variationLatitude = (latLim(2) - latLim(1))/(nLatitudSlices - 1) ;
variationLongitude = (lonLim(2) - lonLim(1))/(nLongitudSlices - 1) ;

regionsLatitude = transpose(latLim(1):variationLatitude:latLim(2)) ;
regionsLongitude = transpose(lonLim(1):variationLongitude:lonLim(2)) ;

[nRowsRegionsLatitude, nColumnsRegionsLatitude] = size(regionsLatitude) ;
[nRowsRegionsLongitude, nColumnsRegionsLongitude] = size(regionsLongitude) ;

nRowsRegionsCoordinates = nRowsRegionsLatitude*nRowsRegionsLongitude ;
nColumnsRegionsCoordinates = nColumnsRegionsLatitude + nColumnsRegionsLongitude ;

regionsCoordinates = zeros(nRowsRegionsCoordinates, nColumnsRegionsCoordinates) ;

k = 0 ;
for i=1:nRowsRegionsLatitude
    for j = 1:nRowsRegionsLongitude
        k = k + 1 ;
        regionsCoordinates(k,:) = [regionsLatitude(i,1), regionsLongitude(j,1)] ;
        
    end
end



%% OBATAINING A VECTOR WITH THE TOTAL EMISSIONS OF CO2 OF EACH REGION

% To get that, a new matrix called "auxiliarMatrix" is going to be used. This matrix 
% is going to report if the dataEdgar point i is inside the region j, and
% in the case that the answer was a yes, the emitted CO2 value is going to
% be save in the postion i,j of this matrix. Using this mechanism, the total CO2 emitted 
% in each region can be got by adding up all the rows of the "auxiliarMatrix" for each one
% of the regions j (regionsEmission vector).

nRowsDataEdgar = gridCoordinatesRows ;

auxiliarMatrix = zeros(nRowsDataEdgar,nRowsRegionsCoordinates) ;

for i=1:nRowsDataEdgar
    for j=1:nRowsRegionsCoordinates
        if dataEdgar(i,1) >= regionsCoordinates(j,1) - variationLatitude/2 && dataEdgar(i,1) < regionsCoordinates(j,1) + variationLatitude/2 && dataEdgar(i,2) >= regionsCoordinates(j,2) - variationLongitude/2 && dataEdgar(i,2) < regionsCoordinates(j,2) + variationLongitude/2
            auxiliarMatrix(i,j) = dataEdgar(i,3) ;

        else
            auxiliarMatrix(i,j) = 0 ;
        end
    end
end

regionsEmission = transpose(sum(auxiliarMatrix,1)) ;
spainEmissions = sum(regionsEmission) ;   

clear auxiliarMatrix

%% OBTAINING SEQUESTRATION DATA

[sequestrationCO2Capacity,depthAquifers] = ExtractingStorageData(regionsCoordinates,nRowsRegionsCoordinates,variationLatitude,variationLongitude) ;

%% IDENTIFYAN POSITION OF EACH REGION (COAST, INLAND, MARITIME)
outputFormat = 'logical' ;
[inlandRegionsMatrix,coastRegionsMatrix,maritimeRegionsMatrix,europeRegions,EU, EUbuf] = GettingLocationRegions(regionsCoordinates,variationLatitude,variationLongitude,latLim,lonLim,outputFormat) ;

% >>> IMPORTANT: always convert to logical mask for indexing
mask = logical(europeRegions);

% Trim ALL arrays to Europe subset
regionsCoordinates       = regionsCoordinates(mask,:);
regionsEmission          = regionsEmission(mask);
sequestrationCO2Capacity = sequestrationCO2Capacity(mask);
depthAquifers            = depthAquifers(mask);

% Trim the classification vectors too (so sizes stay consistent)
inlandRegionsMatrix      = inlandRegionsMatrix(mask);
coastRegionsMatrix       = coastRegionsMatrix(mask);
maritimeRegionsMatrix    = maritimeRegionsMatrix(mask);

%% BUILDING CONNECTED REGIONS SET

nRowsRegionsCoordinates = size(regionsCoordinates,1);


% CALCULATING THE DISTANCE BETWEEN REGIONS

regionsDistance = zeros(nRowsRegionsCoordinates,nRowsRegionsCoordinates) ;
for i=1:nRowsRegionsCoordinates
    for j=1:nRowsRegionsCoordinates
        regionsDistance(i,j) = distance(regionsCoordinates(i,:),regionsCoordinates(j,:),wgs84) /1000 ; % km
    end
end

[connectedRegionsMatrix] = ConnectedRegions(nRowsRegionsCoordinates,regionsDistance) ;

%% BUILDING SURROUNDING REGIONS SET


[surroundingRegionsMatrix] = surroundingRegionsScript(regionsCoordinates,nRowsRegionsCoordinates,variationLatitude,variationLongitude) ;
%% OBTAINING THE SIZE OF EEACH REGION (CELL_SIZE)

intraCellDist = CellSizeCalculation(regionsCoordinates,nRowsRegionsCoordinates,latLim,lonLim,variationLatitude,variationLongitude) ;

%% BUILDING GDX FILE

%egionsAndDataGDX_CR(regionsEmission,regionsCoordinates,sequestrationCO2Capacity,connectedRegionsMatrix,inlandRegionsMatrix,coastRegionsMatrix,maritimeRegionsMatrix,intraCellDist,regionsDistance,depthAquifers,surroundingRegionsMatrix) 
createGDX(regionsEmission, regionsCoordinates, sequestrationCO2Capacity, ...
          connectedRegionsMatrix, inlandRegionsMatrix, coastRegionsMatrix, ...
          maritimeRegionsMatrix, intraCellDist, regionsDistance, depthAquifers);


%% PLOTTING THE REGIONS AND THER RESPECTIVE EMISSIONS

nmeshslices = 100;
xlin = linspace(min(regionsCoordinates(:,2)), max(regionsCoordinates(:,2)),nmeshslices) ;
ylin = linspace(min(regionsCoordinates(:,1)), max(regionsCoordinates(:,1)),nmeshslices) ;

[X,Y] = meshgrid(xlin,ylin) ;

Z = griddata(regionsCoordinates(:,2),regionsCoordinates(:,1),regionsEmission,X,Y,'cubic') ;

figure
mesh(X,Y,Z,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','interp','EdgeColor','interp')
ylim(latLim)
xlim(lonLim) 
ylabel('Lat, °')
xlabel('Lon, °')
zlabel('CO2 emissions, ton/yr')
hold on
plot3(regionsCoordinates(:,2),regionsCoordinates(:,1),regionsEmission,'.','MarkerSize',8,'Color','r')


hold on

land = readgeotable("landareas.shp") ;

geoshow(land,"FaceColor",[1 1 1])

view(2) % 2D view
emissionsColorBar = colorbar;
emissionsColorBar.Label.String = 'CO2 emissions (ton/year)';

Z2 = griddata(regionsCoordinates(:,2),regionsCoordinates(:,1),sequestrationCO2Capacity,X,Y,'cubic') ;

figure
mesh(X,Y,Z2,'EdgeAlpha',0,'FaceAlpha',0.5,'FaceColor','interp','EdgeColor','interp')
ylim(latLim)
xlim(lonLim) 
ylabel('Lat, °')
xlabel('Lon, °')
zlabel('CO2 sequestration capacity, ton CO2')

hold on
plot3(regionsCoordinates(:,2),regionsCoordinates(:,1),sequestrationCO2Capacity+1,'.','MarkerSize',8,'Color','r')

hold on

load coastlines.mat
geoshow(coastlat,coastlon,'Color','k')

view(2) % 2D view
storageCapacityColorBar = colorbar;
storageCapacityColorBar.Label.String = 'Storage Capacity (ton)';


toc

