function   GNSSRprocessing(this,ns,radiusOfEarth)
%% this function:
% -reads location data of a GNSS-R signal receiving CubeSat, or several
% -either computes the locattion of a GNSS-R sending satellite using a Kepler
% propagator, or reads location data of a GNSS-R sending satellite, or several
% -computes the location of the specular point
% -does some basic visualization in 2D or 3D
% does some basic statistical analysis of the specular point location data
%% input variables
% ns: number of CubeSat satellites [-]
% radiusOf Earth [m]
%% output variables
% none

%% begin input section----------------------------------------------------------
%% what visualization and statiscal analysis shall be done? 1=yes, 0=no
plotSPlocationIn3D=1; 
plotSPlocationIn2D=0; 
plotStats=1;          

%% do you want to use the built-in Kepler propagator for the GNSS satellites? 1=yes, 0=no.
%% if you don't, you need to provide the position of the GNSS satellite as input files
useKeplerProp4GNSS=1;

%% define here the size of the specular point. Here, they are only visualized as circles.
%% typically they are ellipses though. For more information see for instance: 
%% N. Rodriguez-Alvarez, E. Podest, K. Jensen, and K. C. McDonald, “Classifying
%% inundation in a tropical wetlands complex with GNSS-R,” Remote Sens., vol. 11, no. 9, 2019.
sizeOfSpecularPoint=10; %% radius of specular point, [km]

%% choose constellation(s): %% 1:Galileo, 2: GPS
%constellations=[1];  %% only Galileo
%constellations=[2];  %% only GPS
%constellations=[1 2];  %% both Galileo and GPS
%constellations=[1 2 3];  %% both Galileo, GPS, Glonass-M
constellations=[1 2 3 4];  %% both Galileo, GPS, Glonass-M, Beidou-2 MEO
%constellations=[1 2 3 4 5];  %% both Galileo, GPS, Glonass-M, Beidou-2 MEO, Beidou-2 non-MEO

%% end input section------------------------------------------------------------


fprintf('\nGNSS processing...');
GNSScpuStartTime = posixtime(datetime('now')); % Posixtime [seconds].

%% read CubeSat position data from file. for data format, see example files
for i=1:ns+1 
  satLLR=readmatrix(strcat('sat',num2str(i-1),'_LLR.csv'));
  timeCubeSat(i,:) = satLLR(:,1);
  latCubeSat(i,:)  = satLLR(:,2);
  lonCubeSat(i,:)  = satLLR(:,3);
  radCubeSat(i,:)  = satLLR(:,4);
end

%% define the GNSS constellation: either propoage using the kepler function or read data file
if useKeplerProp4GNSS %% option 1: compute GNSS location using a kepler function
  %% define initial conditions GNSS constellation
  noOfGNSSsats=0;
  noOfGNSSsatsArray=[];
  inclinationGNSS=[];
  RAANGNSS=[];
  v0GNSS=[];
  altitudeGNSS=[];
  for i=1:size(constellations,2)
    [noOfGNSSsatsPC,inclinationGNSSPC,RAANGNSSPC,v0GNSSPC,altitudeGNSSPC,varargin] = GNSSConstellation(constellations(i));
    noOfGNSSsatsArray=[noOfGNSSsatsArray; noOfGNSSsatsPC];
    noOfGNSSsats=noOfGNSSsats+noOfGNSSsatsPC;
    inclinationGNSS=[inclinationGNSS; inclinationGNSSPC];
    RAANGNSS=[RAANGNSS; RAANGNSSPC];
    v0GNSS=[v0GNSS; v0GNSSPC];
    altitudeGNSS=[altitudeGNSS; altitudeGNSSPC];
    vararginArray(i,1:3)=varargin;
  end  
  
  %% propagate the postion of GNSS-R satellite constellation
  for i=1:size(constellations,2)
    for j=1:noOfGNSSsatsArray(i)
      [timeGNSS(i,j,:),latGNSS(i,j,:),lonGNSS(i,j,:),radGNSS(i,j,:)]=this.keplerPropagation(timeCubeSat(1,:)',timeCubeSat(1,2)-timeCubeSat(1,1),inclinationGNSS(i,j),RAANGNSS(i,j),v0GNSS(i,j),altitudeGNSS(i,j),radiusOfEarth);
      altGNSS(i,j,:)=radGNSS(i,j,:)-radiusOfEarth/1000;
    end
  end
else %% option 2: read GNSS location from data file, this option is not tested
%  for i=1:noOfGNSSsats
%   GNSSsatLLR=readmatrix(strcat('GNSSsat',num2str(i-1),'_LLR.csv'));
%   timeGNSS(i,:) = GNSSsatLLR(:,1);
%   latGNSS(i,:)  = GNSSsatLLR(:,2);
%   lonGNSS(i,:)  = GNSSsatLLR(:,3);
%   altGNSS(i,:)  = GNSSsatLLR(:,4)-radiusOfEarth;  
%  end
  ;
end

fprintf('\nsetting up GNSS constellation time: %s seconds.',num2str(posixtime(datetime('now')) - GNSScpuStartTime));

fprintf('\ncomputing specular point location');

%% compute SP location per each cubesat, per each GNSS constellation and per each GNSS satellite
latSP=zeros(ns+1,size(constellations,2),max(noOfGNSSsatsArray),size(timeCubeSat,2));
lonSP=zeros(ns+1,size(constellations,2),max(noOfGNSSsatsArray),size(timeCubeSat,2));
for i=2:ns+1
  for j=1:size(constellations,2)
    for k=1:noOfGNSSsatsArray(j)
      [latSP(i,j,k,:), lonSP(i,j,k,:)]  =  computeSPlocation(timeCubeSat(i,:),latCubeSat(i,:)',lonCubeSat(i,:)',radCubeSat(i,:)',squeeze(latGNSS(j,k,:))',squeeze(lonGNSS(j,k,:))',squeeze(altGNSS(j,k,:))'+radiusOfEarth/1000,radiusOfEarth/1000);
      if i==2 && j==1 && k==1%% screenoutput on progress
        fprintf('\ncomputing specular points ')
        fprintf('%4.0f/%4.0f',(i-2)*size(constellations,2)*max(noOfGNSSsatsArray) +(j-1)*max(noOfGNSSsatsArray)+k,ns*size(constellations,2)*max(noOfGNSSsatsArray));
      else
        %fprintf('\b\b\b\b\b\b\b\b\b');
        fprintf('\n%4.0f/%4.0f',(i-2)*size(constellations,2)*max(noOfGNSSsatsArray) +(j-1)*max(noOfGNSSsatsArray)+k,ns*size(constellations,2)*max(noOfGNSSsatsArray));
      end
    end
  end
end

save('latSP.mat','latSP');
save('lonSP.mat','lonSP');
fprintf('\ncomputing specular point location time: %s seconds.',num2str(posixtime(datetime('now')) - GNSScpuStartTime));

%% 3D plot: position of the CubeSat, the GNSS satellites and the specular point
if plotSPlocationIn3D 
  fprintf('\ndisplaying in 3D');
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
  %% display position of the CubeSats
  for i=2:ns+1
    plotm(latCubeSat(i,:),lonCubeSat(i,:),radCubeSat(i,:)-radiusOfEarth/1000,'Color',[0 0 0]);hold on;
  end
  %% display position of GNSS satellites
  for i=1:size(constellations,2)
    for j=1:noOfGNSSsatsArray(i)
      plotm(squeeze(latGNSS(i,j,:)),squeeze(lonGNSS(i,j,:)),squeeze(altGNSS(i,j,:)),'Color', vararginArray{i,2});hold on;
    end
  end
  %% display location of SP
  for i=2:ns+1
    for j=1:size(constellations,2)
      for k=1:noOfGNSSsatsArray(j)
        p1=plotm(squeeze(latSP(i,j,k,:))',squeeze(lonSP(i,j,k,:))',ones(1,size(lonSP,4)),'.','Color', vararginArray{j,2});hold on;
        %circlem(squeeze(latSP(i,j,:)),squeeze(lonSP(i,j,:)),10000*ones(size(lonSP,3),1));hold on; %% does not work at the moment
      end
      plotname(j)=p1;
      legendname(j)=vararginArray(j,3);
    end
  end
  legend(plotname,legendname);
  view(90,0)
  fprintf('\ndisplaying 3D time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
end %% plot3D

%%2D plot: position of the CubeSat and the specular point
if plotSPlocationIn2D 
  fprintf('\ndisplaying in 2D');
  figure
  ax = worldmap('World');
  setm(ax, 'Origin', [0 0 0])
  land = shaperead('landareas', 'UseGeoCoords', true);
  geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
  lakes = shaperead('worldlakes', 'UseGeoCoords', true);
  geoshow(lakes, 'FaceColor', 'blue')
  rivers = shaperead('worldrivers', 'UseGeoCoords', true);
  geoshow(rivers, 'Color', 'blue')
  for i=2:ns+1
    %% plot subsatellite points
    plotm(latCubeSat(i,:),lonCubeSat(i,:),'.k');hold on;
    %% plot specular point for each GNSS constellation
    %for j=1:noOfGNSSsats %% define color of specular point for each constellation
    for j=1:size(constellations,2) %% define color of specular point for each constellation
      for k=1:noOfGNSSsatsArray(j)
        p1=circlem(squeeze(latSP(i,j,k,:))',squeeze(lonSP(i,j,k,:))',sizeOfSpecularPoint*ones(size(lonSP,4),1),{'EdgeColor',vararginArray{j,2}});hold on;
      end
      if i==2 && j==1
        fprintf('\ndisplaying specular points ')
        fprintf('%4.0f/%4.0f',(i-2)*size(constellations,2)+j,ns*size(constellations,2));
      else
        fprintf('\b\b\b\b\b\b\b\b\b');
        fprintf('%4.0f/%4.0f',(i-2)*size(constellations,2)+j,ns*size(constellations,2));
      end
      plotname(j)=p1;
      legendname(j)=vararginArray(j,3);
    end
  end
  legend(plotname,legendname);
  fprintf('\ndisplaying 2D time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
end %% plot 2D

%% basic statistical analysis and visualization of the location of the specular point
if plotStats
 % totalLat=zeros(size(latSP,1)*size(latSP,3)*size(latSP,3),1);
 % totalLon=zeros(size(latSP,1)*size(latSP,3)*size(latSP,3),1);
 m=1;
 for i=2:size(latSP,1)
   for j=1:size(latSP,2)
     for k=1:size(latSP,3)
       for l=1:size(latSP,4)
         totalLat(m)=latSP(i,j,k,l);
         totalLon(m)=wrapTo180(lonSP(i,j,k,l));
         m=m+1;
       end
     end
   end
 end
 
 figure
 latGranularity=1; %% [deg]
 lonGranularity=1; %% [deg]
 hist3([totalLon' totalLat'],'Ctrs',{-180:lonGranularity:180 -90:latGranularity:90})
end %% plot statistics


%% selection of sp, e.g.: 1) one point per geosurface (world, land, europe), compute results/effort 2) compute EO campaign optimcally distributed over formation members


fprintf('\n Total GNSS computation time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
end %% function GNSS-R processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
function [noOfGNSSsats,inclinationGNSS,RAANGNSS,trueAnomalyGNSS,altitudeGNSS,varargin]=GNSSConstellation(constellation)
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
  case 1 %% Galileo, source: wikipedia, not accounting for spares
    noOfPlanes=3;
    noOfGNSSsatsPerPlane=8;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=56*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=23.222e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=77.632;                        %% [deg]
    trueAnomalyShift=0;                      %% [deg]
    varargin={'constellationColor',[1 0 0],'Galileo'}; %red
  case 2 %% GPS, source: wikipedia, not accounting for spares
    noOfPlanes=6;
    noOfGNSSsatsPerPlane=4;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=55*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=20.2e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=46.7112;                       %% [deg]
    trueAnomalyShift=0;                      %% [deg]
    varargin={'constellationColor',[1 1 0],'GPS'}; %yellow
  case 3 %% Glonass-M, source: wikipedia, not accounting for spares
    noOfPlanes=3;
    noOfGNSSsatsPerPlane=8;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=64.8*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=19.1e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=46.7112;                       %% [deg]
    trueAnomalyShift=0;                      %% [deg]
     varargin={'constellationColor',[0 0 1],'Glonass-M'}; %blue
  case 4 %% Beidou-2, only MEO satellites, source: https://gssc.esa.int/navipedia/index.php/BeiDou_Space_Segment, not accounting for spares
    noOfPlanes=3;
    noOfGNSSsatsPerPlane=8;
    noOfGNSSsats=noOfGNSSsatsPerPlane*noOfPlanes;
    inclinationGNSS=55*ones(1,noOfGNSSsats); %% [deg]
    altitudeGNSS=27.878e6*ones(1,noOfGNSSsats);  %% [km]
    RAANshift=46.7112;                       %% [deg]
    trueAnomalyShift=0;                      %% [deg]
    varargin={'constellationColor',[0 1 0],'Beidou MEO'};%green 
  case 5 %% Beidou-2, non-MEO satellites
    %TBD
    ;
    varargin={'constellationColor',[0.5 0.5 0],'Beidou non-MEO'};%dirty green
  case 6 %% Indian
    %TBD
    varargin={'constellationColor',[0 1 1],'Indian GNSS'}; %light blue
  case 7 %% Japanese
    %TBD
    varargin={'constellationColor',[1 0 1],'Japanese GNSS'}; %purple
  otherwise
    fprintf('\nerror setting up GNSS constellation');
    varargin={'constellationColor',[0 0.5 0.5]}; %	turquoise
end

RAANGNSS=zeros(1,noOfPlanes*noOfGNSSsatsPerPlane);  %% right ascension of ascending node [deg]
trueAnomalyGNSS=zeros(1,noOfPlanes*noOfGNSSsatsPerPlane);    %% true anomaly[deg]

%% equistant distribution of orbital planes and satellites on the planes
for i=1:noOfPlanes
  for j=1:noOfGNSSsatsPerPlane
    RAANGNSS((i-1)*noOfGNSSsatsPerPlane+j)=360/noOfPlanes*(i-1)+RAANshift;
    trueAnomalyGNSS((i-1)*noOfGNSSsatsPerPlane+j)=360/noOfGNSSsatsPerPlane*(j-1)+trueAnomalyShift;
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
%% does this work?:
% https://nl.mathworks.com/help/map/ref/ellipse1.html
% if so implement ellipses according to:
% "Classifying Inundation in a TropicalWetlands Complex with GNSS-R" by Nereida Rodriguez-Alvarez 1,*, Erika Podest 2, Katherine Jensen 3,4 and Kyle C. McDonald 2,3,4
%% Calculate circle coordinates:
[circlelat,circlelon] = scircle1(lat,lon,radius,[],earthRadius(units));


%%  Plot and format circle(s): 
h = patchm(circlelat,circlelon,'k','facecolor','none');
 if nargin>3 && ~isempty(varargin)
     %set(h,varargin{:})
     set(h,varargin{1}{1},varargin{1}{2})
 end
%% Clean up: 
if nargout==0
    clear h
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       