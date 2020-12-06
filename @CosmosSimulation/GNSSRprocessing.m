function   GNSSRprocessing(this,ns,radiusOfEarth)
fprintf('\n GNSSR processing function')

plotSPlocationIn3D=1; %% works only if runSPfinder is enabled
plotSPlocationIn2D=1; %% works only if runSPfinder is enabled
plotStats=1;

sizeOfSpecularPoint=10; %% radius, [km]

%% choose constellation(s): %% 1:Galileo, 2: GPS
%constellations=[1]; 
constellations=[2];
%constellations=[1 2]; 

%% read data from file. for data format, see example files
for i=1:ns+1 
  satLLR=readmatrix(strcat('sat',num2str(i-1),'_LLR.csv'));
  timeCubeSat(i,:) = satLLR(:,1);
  latCubeSat(i,:)  = satLLR(:,2);
  lonCubeSat(i,:)  = satLLR(:,3);
  radCubeSat(i,:)  = satLLR(:,4);
end

fprintf('\nGNSS processing...');
GNSScpuStartTime = posixtime(datetime('now')); % Posixtime [seconds].

%% define initial conditions GNSS constellation
noOfGNSSsats=0;
noOfGNSSsatsArray=[];
inclinationGNSS=[];
RAANGNSS=[];
v0GNSS=[];
altitudeGNSS=[];

for i=1:size(constellations,2)
  [noOfGNSSsatsPC,inclinationGNSSPC,RAANGNSSPC,v0GNSSPC,altitudeGNSSPC] = GNSSConstellation(constellations(i));
  noOfGNSSsatsArray=[noOfGNSSsatsArray; noOfGNSSsatsPC];
  noOfGNSSsats=noOfGNSSsats+noOfGNSSsatsPC;
  inclinationGNSS=[inclinationGNSS; inclinationGNSSPC];
  RAANGNSS=[RAANGNSS; RAANGNSSPC];
  v0GNSS=[v0GNSS; v0GNSSPC];
  altitudeGNSS=[altitudeGNSS; altitudeGNSSPC];
end

%% propagate GNSS constellation
for i=1:noOfGNSSsats
  [timeGNSS(i,:),latGNSS(i,:),lonGNSS(i,:),radGNSS(i,:)]=this.keplerPropagation(timeCubeSat(1,:)',timeCubeSat(1,2)-timeCubeSat(1,1),inclinationGNSS(i),RAANGNSS(i),v0GNSS(i),altitudeGNSS(i),radiusOfEarth);
  altGNSS(i,:)=radGNSS(i,:)-radiusOfEarth/1000;
end
fprintf('\nsetting up GNSS constellation time: %s seconds.',num2str(posixtime(datetime('now')) - GNSScpuStartTime));

%% compute SP location per each cubesat and per each GNSS satellite
latSP=zeros(ns+1,noOfGNSSsats,size(timeCubeSat,2));
lonSP=zeros(ns+1,noOfGNSSsats,size(timeCubeSat,2));
for i=2:ns+1
  for j=1:noOfGNSSsats
    [latSP(i,j,:), lonSP(i,j,:)]=computeSPlocation(timeCubeSat(i,:),latCubeSat(i,:)',lonCubeSat(i,:)',radCubeSat(i,:)',latGNSS(j,:),lonGNSS(j,:),altGNSS(j,:)+radiusOfEarth/1000,radiusOfEarth/1000);
    if i==2 && j==1
      fprintf('\ncomputing specular points ')
      fprintf('%4.0f/%4.0f',(i-2)*noOfGNSSsats+j,ns*noOfGNSSsats);
    else
      fprintf('\b\b\b\b\b\b\b\b\b');
      fprintf('%4.0f/%4.0f',(i-2)*noOfGNSSsats+j,ns*noOfGNSSsats);
    end
  end
end
fprintf('\ncomputing specular point location time: %s seconds.',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
%latSP
%size(latSP)
if plotSPlocationIn3D %% 3D plot
  %% set-up GIS
  grs80 = referenceEllipsoid('grs80','km');
  load topo
  figure('Renderer','opengl')
  ax = axesm('globe','Geoid',grs80,'Grid','off','GLineWidth',1,'GLineStyle','-', 'Gcolor',[0.9 0.9 0.1],'Galtitude',100);
  ax.Position = [0 0 1 1];
  axis equal off
  geoshow(topo,topolegend,'DisplayType','texturemap')
  demcmap(topo)
  land = shaperead('landareas','UseGeoCoords',true);
  plotm([land.Lat],[land.Lon],'Color','black')
  rivers = shaperead('worldrivers','UseGeoCoords',true);
  plotm([rivers.Lat],[rivers.Lon],'Color','blue')
  
  %% display position satellites
  for i=2:ns+1
    plotm(latCubeSat(i,:),lonCubeSat(i,:),radCubeSat(i,:)-radiusOfEarth/1000,'Color',[1 0  0]);hold on;
  end
  %% display position of GNSS satellites
  for i=1:noOfGNSSsats
    plotm(latGNSS(i,:),lonGNSS(i,:),altGNSS(i,:),'Color',[0 1 0]);hold on;
  end
  %% display location of SP
  for i=2:ns+1
    for j=1:noOfGNSSsats
      plotm(squeeze(latSP(i,j,:))',squeeze(lonSP(i,j,:))',ones(1,size(lonSP,3)),'.');hold on;
      %circlem(squeeze(latSP(i,j,:)),squeeze(lonSP(i,j,:)),10000*ones(size(lonSP,3),1));hold on; %% does not work at the moment
    end
  end
  view(90,0)
  fprintf('\n displaying 3D time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
end

if plotSPlocationIn2D %%2D plot
  figure
  ax = worldmap('World');
  setm(ax, 'Origin', [0 0 0])
  land = shaperead('landareas', 'UseGeoCoords', true);
  geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
  lakes = shaperead('worldlakes', 'UseGeoCoords', true);
  geoshow(lakes, 'FaceColor', 'blue')
  rivers = shaperead('worldrivers', 'UseGeoCoords', true);
  geoshow(rivers, 'Color', 'blue')
  SPcolor=rand(noOfGNSSsats,3);
  
  for i=2:ns+1
    %% plot subsatellite points
    plotm(latCubeSat(i,:),lonCubeSat(i,:),'.');hold on;
    for j=1:noOfGNSSsats
      if j<=noOfGNSSsatsArray(1)
        varargin={'EdgeColor',[1 0 0]};
      else
        varargin={'EdgeColor',[0 0 0]};
      end
      %varargin={'EdgeColor',SPcolor(j,:)};
      varargin=[];
      circlem(squeeze(latSP(i,j,:))',squeeze(lonSP(i,j,:))',sizeOfSpecularPoint*ones(size(lonSP,3),1),varargin);hold on;
      if i==2 && j==1
        fprintf('\ndisplaying specular points ')
        fprintf('%4.0f/%4.0f',(i-2)*noOfGNSSsats+j,ns*noOfGNSSsats);
      else
        fprintf('\b\b\b\b\b\b\b\b\b');
        fprintf('%4.0f/%4.0f',(i-2)*noOfGNSSsats+j,ns*noOfGNSSsats);
      end
    end
  end
  fprintf('\n displaying 2D time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
  
end

if plotStats
 % totalLat=zeros(size(latSP,1)*size(latSP,3)*size(latSP,3),1);
 % totalLon=zeros(size(latSP,1)*size(latSP,3)*size(latSP,3),1);
  l=1;
  for i=2:size(latSP,1)
    for j=1:size(latSP,2)
      for k=1:size(latSP,3)
        totalLat(l)=latSP(i,j,k);
        totalLon(l)=lonSP(i,j,k);        
        l=l+1;
      end
    end
  end
figure
hist3([totalLon' totalLat'],'Ctrs',{0:5:360 -90:5:90})
 
end

fprintf('\n Total GNSS computation time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));

end

function [latSP, lonSP]=computeSPlocation(time,latCubeSat,lonCubeSat,radCubeSat,latGNSS,lonGNSS,radGNSS,rC)
%% based on:
%% A New Approach to Determine the Specular Point of Forward Reflected GNSS Signals
%% Benjamin John Southwell Dolby Laboratories, Inc. Andrew G Dempster UNSW Sydney
%%%% input variables:
%% vizgridtime
%% latCubeSat
%% lonCubeSat
%% radCubeSat
%% timeGNSS
%% latGNSS
%% lonGNSS
%% radGNSS
%% radiusOfEarth
%%%% output variables:
%% latSP
%% lonSP

latSP=zeros(1,size(time,1));
lonSP=zeros(1,size(time,1));

for i=1:size(latGNSS,2) %% compute betas for each lat data point
  T=radGNSS(i);
  R=radCubeSat(i);
  %% compute thetaR, thetaT
  if latGNSS(i)> 270 && latCubeSat(i)<90
    THETAlat=latGNSS(i)-latCubeSat(i)-360;
  elseif latGNSS(i) < 90 && latCubeSat(i) > 270
    THETAlat=latGNSS(i)-latCubeSat(i)+360;
  else
    THETAlat=latGNSS(i)-latCubeSat(i);
  end
  
  myfunLat = @(betaLat) asind(rC/T*sind(betaLat))+asind(rC/R*sind(betaLat))+THETAlat+2*betaLat-360;
  funLat = @(betaLat) myfunLat(betaLat);
  betaLat = fzero(funLat,180-THETAlat);
  betaLatOld=betaLat;
  thetaSPLat(i)=180-betaLat-asind(rC/T*sind(betaLat));
  latSPTemp(i)=latGNSS(i)-thetaSPLat(i);
  %fprintf('\n thetaSPLat %f betaLat %f latSPTemp %f: ',thetaSPLat(i),betaLat,latSPTemp(i))
  
  if lonGNSS(i)> 270 && lonCubeSat(i)<90
    THETAlon=lonGNSS(i)-lonCubeSat(i)-360;
  elseif lonGNSS(i) < 90 && lonCubeSat(i) > 270
    THETAlon=lonGNSS(i)-lonCubeSat(i)+360;
  else
    THETAlon=lonGNSS(i)-lonCubeSat(i);
  end
  
  %fprintf('\n lon1 %f %f %f: ',lonGNSSnewGrid(i),lonRec(i),THETAlon)
  myfunLon = @(betaLon) asind(rC/T*sind(betaLon))+asind(rC/R*sind(betaLon))+THETAlon+2*betaLon-360;
  funLon = @(betaLon) myfunLon(betaLon);
  betaLon = fzero(funLon,180-THETAlon);
  betaLonOld=betaLon;
  thetaSPLon(i)=180-betaLon-asind(rC/T*sind(betaLon));
  lonSPTemp(i)=lonGNSS(i)-thetaSPLon(i);
  %fprintf('\n lon2 %f %f %f : ',thetaSPLon(i),betaLon,lonSPTemp(i))
  
  if isfinite(latSPTemp(i)) && isfinite(lonSPTemp(i))
    if abs(latSPTemp(i))>90
      fprintf('\n error in latitude: latSPTemp(i) %f ',latSPTemp(i));
    else
      latSP(i)=latSPTemp(i);
    end
    lonSP(i)=wrapTo360(lonSPTemp(i));
  else
    latSP(i)=0;
    lonSP(i)=0;
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    

function [noOfGNSSsats,inclinationGNSS,RAANGNSS,trueAnomalyGNSS,altitudeGNSS]=GNSSConstellation(constellation)
%% this function defines a constellation
%% input variables
% n/a
%% output variables
% noOfGNSSsats
% inclinationGNSS
% RAANGNSS
% v0GNSS
% altitudeGNSS
%% internal variables
% noOfPlanes
% noOfGNSSsatsPerPlane
% noOfGNSSsats
% RAANGNSS
% v0GNSS


switch constellation
  case 1 %% Galileo
    noOfPlanes=3;
    noOfGNSSsatsPerPlane=10;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=55*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=20e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=77.632;                        %% [deg]
    trueAnomalyShift=0;                      %% [deg]
  case 2 %% GPS
    noOfPlanes=6;
    noOfGNSSsatsPerPlane=5;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=55*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=20e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=46.7112;                       %% [deg]
    trueAnomalyShift=0;                      %% [deg]
  otherwise
    fprintf('\nerror setting up GNSS constellation');
end

RAANGNSS=zeros(noOfPlanes*noOfGNSSsatsPerPlane);  %% right ascension of ascending node [deg]
trueAnomalyGNSS=zeros(noOfPlanes*noOfGNSSsatsPerPlane);    %% true anomaly[deg]

%% equistant distribution of orbital planes and satellites
for i=1:noOfPlanes
  for j=1:noOfGNSSsatsPerPlane
    RAANGNSS((i-1)*j+j)=360/noOfPlanes*(i-1)+RAANshift;
    trueAnomalyGNSS((i-1)*j+j)=360/noOfGNSSsatsPerPlane*(j-1)+trueAnomalyShift;
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
         

function [h,circlelat,circlelon] = circlem(lat,lon,radius,varargin)
%CIRCLEM draws circles on maps. 
% 
%% Syntax
% 
%  circlem(lat,lon,radius)
%  circlem(...,'units',LengthUnit)
%  circlem(...,'PatchProperty',PatchValue)
%  h = circlem(...)
%  [h,circlelat,circlelon] = circlem(...)
% 
%% Description 
%
% circlem(lat,lon,radius) draws a circle or circles of radius or radii
% given by radius centered at lat, lon, where radius, lat, and
% lon may be any combination of scalars, vectors, or MxN array. All non-
% scalar inputs must have matching dimensions. 
% 
% circlem(...,'units',LengthUnit) specifies a length unit of input
% radius. See validateLengthUnit for valid units. Default unit is
% kilometers. 
% 
% circlem(...,'PatchProperty',PatchValue) specifies patch properties such
% as edgecolor, facecolor, facealpha, linewidth, etc. 
% 
% h = circlem(...) returns the patch handle of plotted circle(s). 
%  
% [h,circlelat,circlelon] = circlem(...) also returns arrays of latitudes 
% and longitudes corresponding to the outline of each circle drawn.  Each
% "circle" is in actuality a polygon made of 100 lat/lon pairs.
% 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas Institute 
% for Geophysics (UTIG) on October 14, 2014.  
%
% See also scircle1, distance, validateLengthUnit
%% Initial input checks: 
assert(license('test','map_toolbox')==1,'circlem requires Matlab''s Mapping Toolbox.')    
try
    MapIsCurrent = ismap; 
    assert(MapIsCurrent==1,'The circlem function requires you to initialize a map first.') 
catch
    error('The circlem function requires you to initialize a map first.') 
end    
assert(nargin>2,'circlem function requires at least three inputs--latitude(s), longitude(s), and radius(ii)')
assert(isnumeric(lat)==1,'Input latitude(s) must be numeric.');
assert(isnumeric(lon)==1,'Input longitude(s) must be numeric.');
assert(max(abs(lat))<=90,'Latitudes cannot exceed 90 degrees North or South.') 
assert(max(abs(lon))<=360,'Longitudes cannot exceed 360 degrees North or South.') 
assert(isnumeric(radius)==1,'Radius must be numeric.') 
%% Declare units
units = 'km'; % kilometers by default
tmp = strncmpi(varargin,'unit',4); 
if any(tmp)
    units = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
end
%% Reshape inputs as needed
% columnate:
lat = lat(:); 
lon = lon(:); 
radius = radius(:); 
% How many circles do we need to make? 
NumCircles = max([length(lat) length(lon) length(radius)]); 
% Make vectors all the same lengths so scircle1 will be happy:  
if length(lat)<NumCircles
    assert(isscalar(lat)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    lat = lat*ones(NumCircles,1); 
end
if length(lon)<NumCircles
    assert(isscalar(lon)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    lon = lon*ones(NumCircles,1); 
end
if length(radius)<NumCircles
    assert(isscalar(radius)==1,'It seems that the inputs to circlem have too many different sizes. Lat, lon, and radius can be any combination of scalars, vectors, or 2D grids, but all nonscalar inputs must be the same size.')
    radius = radius*ones(NumCircles,1); 
end
%% Calculate circle coordinates:
[circlelat,circlelon] = scircle1(lat,lon,radius,[],earthRadius(units));

%% does this work?:
% https://nl.mathworks.com/help/map/ref/ellipse1.html
% if so implement ellipses according to:
% "Classifying Inundation in a TropicalWetlands Complex with GNSS-R" by Nereida Rodriguez-Alvarez 1,*, Erika Podest 2, Katherine Jensen 3,4 and Kyle C. McDonald 2,3,4

%%  Plot and format circle(s): 
varargin=[];
h = patchm(circlelat,circlelon,'k','facecolor','none');
 if nargin>3 && ~isempty(varargin)
     set(h,varargin{:})
 end
%% Clean up: 
if nargout==0
    clear h
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

