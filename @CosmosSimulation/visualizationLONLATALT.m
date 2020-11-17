function visualizationLONLATALT(this, vizScale, ns, altitude)
%% input:
% ns: number of satellites
% altitude: altitude used for visualization. The actual altitude may change for long periods of
% operations. However, this will be so small that it is irrelevant for visualization
%% output:
% nonef
% however, files are written that are used by cosmosVIZ

% Get path to telemetry files from CosmosSimulation object.
tmFolderPath = this.TelemetryPath;

%% switches
plotLatLonIn2D=0;
writeLLRRPYData=0;
runSPfinder=1;
plotSPlocationIn3D=0;
plotSPlocationIn2D=1;

%% global parameter for formation
inclination=0; %% [-90, 90] [deg]  %% equatorial
%inclination=35; %% [-90, 90] [deg] %% CYGNSS
%inclination=52; %% [-90, 90] [deg] %% ISS
%inclination=90; %% [-90, 90] [deg] %% polar
%inclination=97; %% [-90, 90] [deg] %% (approximately) SSO
radiusOfEarth=6371000;          %% [m]
RAAN=0; %%RAAN    = input(' Right Ascension of Ascendent Node    [  0,360[    RAAN   [deg] = ');
v0=0;   %%v0      = input(' True anomaly at the departure        [  0,360[    v0     [deg] = ');

% Set path for LLR-RPY files.
rpyFolderPath = 'LLR-RPY';
[~, ~, ~] = mkdir(rpyFolderPath); % [status, msg, msgID]
delete(strcat(rpyFolderPath,filesep,'*')); % delete old files in folder.

% improved the GNSS vizualization. in particular now the sp size is vizualized in 2D
keplerStepSize=20; %% [s]
sizeOfSpecularPoint=10; % radius, [km]
constellations=[1 2]; %% 1:Galileo, 2: GPS


%% read data from telemery files
% JT: this works only if the telemetry data of all satellites is equal in size 
% and sync'ed. Sooner or later, this needs to become more versatile allowing individual times for each satellite
for i=1:ns
  fileControlVectorTM = strcat(tmFolderPath,filesep,'ControlVectorTM',num2str(i),'.csv');
  fileForceVectorTM = strcat(tmFolderPath,filesep,'ForceVectorTM',num2str(i),'.csv');
  fileSatPositionTM = strcat(tmFolderPath,filesep,'SatPositionTM',num2str(i),'.csv');
  fileTimeVectorTM = strcat(tmFolderPath,filesep,'TimeVectorTM',num2str(i),'.csv');
  fileSatStatesTM = strcat(tmFolderPath,filesep,'SatStatesTM',num2str(i),'.csv');
  
  tempTime=readmatrix(fileTimeVectorTM);  
  tempSatStates=readmatrix(fileSatStatesTM); 
  if i==1
      timeSteps=size(tempTime,1);
      cosmosTime=zeros(timeSteps,ns);
      sstXViz =zeros(ns,timeSteps);
      sstYViz =zeros(ns,timeSteps);
      sstZViz =zeros(ns,timeSteps);
      roll =zeros(ns,timeSteps);
      pitch=zeros(ns,timeSteps);
      yaw  =zeros(ns,timeSteps);
  end
  cosmosTime(:,i)=tempTime(:);
  sstX(i,:)  =tempSatStates(:,1)';
  sstY(i,:)  =tempSatStates(:,2)';
  sstZ(i,:)  =tempSatStates(:,3)';
  roll(i,:)  =tempSatStates(:,7)';
  pitch(i,:) =tempSatStates(:,8)';
  yaw(i,:)   =tempSatStates(:,9)';
end


%% rotate relative coordinate system to align x+ with flight direction (inclination)
for i=1:ns
  for j=1:size(sstXViz,2)
    sstTemp=rotz(inclination-90)*[sstX(i,j) sstY(i,j) 0]';
  end
  sstX(i,j)=sstTemp(1);
  sstY(i,j)=sstTemp(2);
end
  
%% compute orbit
[vizTime,latitude,longitude,radius] = kepler(cosmosTime,keplerStepSize,inclination,RAAN,v0,altitude,radiusOfEarth);

%% interpolate relative position on visualization time steps
for i=1:ns
  %% interpolate on viz-time-grid
  sstXvizTime(i,:)      =interp1(cosmosTime(:,i),sstX(i,:),vizTime);
  sstYvizTime(i,:)      =interp1(cosmosTime(:,i),sstY(i,:),vizTime);
  sstZvizTime(i,:)      =interp1(cosmosTime(:,i),sstZ(i,:),vizTime);
  
  %% interpolate on viz-time-grid and apply visualization scaling
  sstXvizTimeScaled(i,:)=interp1(cosmosTime(:,i),sstX(i,:)*vizScale,vizTime);
  sstYvizTimescaled(i,:)=interp1(cosmosTime(:,i),sstY(i,:)*vizScale,vizTime);
  sstZvizTimeScaled(i,:)=interp1(cosmosTime(:,i),sstZ(i,:)*vizScale,vizTime);
  
  %% interpolate Euler angles
  rollVizTime(i,:)      =interp1(cosmosTime(:,i),squeeze(roll(i,:)),vizTime);
  pitchVizTime(i,:)     =interp1(cosmosTime(:,i),squeeze(pitch(i,:)),vizTime);
  yawVizTime(i,:)       =interp1(cosmosTime(:,i),squeeze(yaw(i,:)),vizTime);
end
      
%% formation centerpoint location
lat(1,:)      = latitude;
lon(1,:)      = longitude;
rad(1,:)      = radius;
%% formation centerpoint location for scaled variables
latScaled(1,:)      = latitude;
lonScaled(1,:)      = longitude;
radScaled(1,:)      = radius;

%% formation centerpoint
rollVizTime         = [zeros(1,size(pitchVizTime,2)); rollVizTime];
pitchVizTime        = [zeros(1,size(pitchVizTime,2)); pitchVizTime];
yawVizTime          = [zeros(1,size(pitchVizTime,2)); yawVizTime];

%% off set of the formation satellites
for i=1:size(vizTime,2)
  for j=1:ns
    %% scaled for video visualization
    localTransformedCoSystemScaled = rotz(inclination) * [sstXvizTimeScaled(j,i) ; sstYvizTimescaled(j,i); sstZvizTimeScaled(j,i)];
    LonOffsetAngleScaled(j,i)      = asind( localTransformedCoSystemScaled(1)/1000 / (radScaled(1,i) + localTransformedCoSystemScaled(3)/1000) );
    LatOffsetAngleScaled(j,i)      = asind( localTransformedCoSystemScaled(2)/1000 / (radScaled(1,i) + localTransformedCoSystemScaled(3)/1000) );
    %% not scaled for 2D and 3 visualization
    localTransformedCoSystem       = rotz(inclination) * [sstXvizTime(j,i) ; sstYvizTime(j,i); sstZvizTime(j,i)];
    LonOffsetAngle(j,i)            = asind( localTransformedCoSystem(1)/1000 / (rad(1,i) + localTransformedCoSystem(3)/1000) );
    LatOffsetAngle(j,i)            = asind( localTransformedCoSystem(2)/1000 / (rad(1,i) + localTransformedCoSystem(3)/1000) );
  end
end
  
%% Lat, Lon, Rad of formation satellites
for i=1:size(vizTime,2)
  for j=2:ns+1
    
    lonScaled(j,i)     = wrapTo360(lonScaled(1,i)+LonOffsetAngleScaled(j-1,i));         % Longitude            [deg]
    latScaled(j,i)     = latScaled(1,i)+LatOffsetAngleScaled(j-1,i);                     % Latitude             [deg]
    radScaled(j,i)     = radScaled(1,i)'+sstZvizTimeScaled(j-1,i)/1000;                      % radius               [km]
    
    lon(j,i)           = wrapTo360(lon(1,i)+LonOffsetAngle(j-1,i));                     % Longitude            [deg]
    lat(j,i)           = lat(1,i)+LatOffsetAngle(j-1,i);                                 % Latitude             [deg]
    rad(j,i)           = rad(1,i)'+sstZvizTime(j-1,i)/1000;                                  % radius               [km]

      %{
      if i==1 %% 1-point, northward
          Lat(j,i)     = Lat(1,i)+ramOffsetAngle(j-1,i)/pi*180;              % Latitude             [deg]
          Lon(j,i)     = wrapTo360(Lon(1,i)+planeOffsetAngle(j-1,i)/pi*180); % Longitude            [deg]
          Rad(j,i)     = Rad(i)'+sstzvizgrid(j-1,i)/1000;                    % radius                [km]
        elseif Lat(1,i+1)>Lat(1,i) %% northward
          Lat(j,i)     = Lat(1,i)+ramOffsetAngle(j-1,i)/pi*180;              % Latitude             [deg]
          Lon(j,i)     = wrapTo360(Lon(1,i)+planeOffsetAngle(j-1,i)/pi*180); % Longitude            [deg]
          Rad(j,i)     = Rad(i)'+sstzvizgrid(j-1,i)/1000;                    % radius                [km]
        elseif Lat(1,i+1)<Lat(1,i) %% southward
          Lat(j,i)     = Lat(1,i)-ramOffsetAngle(j-1,i)/pi*180;      % Latitude             [deg]
          Lon(j,i)     = wrapTo360(Lon(1,i)+planeOffsetAngle(j-1,i)/pi*180);    % Longitude            [deg]
          Rad(j,i)     = Rad(i)'+sstzvizgrid(j-1,i)/1000;       % radius                [km]
          pitchvizgrid(j-1,i)= pitchvizgrid(j-1,i)-180;
        elseif i==size(vizgridtime,2) %% last point
            if Lat(1,i)>Lat(1,i-1) %% northward
              Lat(j,i)    = Lat(1,i)+ramOffsetAngle(j-1,i)/pi*180;      % Latitude             [deg]
              Lon(j,i)    = wrapTo360(Lon(1,i)+planeOffsetAngle(j-1,i)/pi*180);    % Longitude            [deg]
              Rad(j,i)    = Rad(i)'+sstzvizgrid(j-1,i)/1000;       % radius                [km]
            elseif Lat(1,i)<Lat(1,i-1) %% southward
              Lat(j,i)     = Lat(1,i)-ramOffsetAngle(j-1,i)/pi*180;      % Latitude             [deg]
              Lon(j,i)     = wrapTo360(Lon(1,i)-planeOffsetAngle(j-1,i)/pi*180);    % Longitude            [deg]
              Rad(j,i)     = Rad(i)'+sstzvizgrid(j-1,i)/1000;       % radius                [km]
            end
        else
          %send(dq,'error visualizationLONLATALT');
        end
      
      %}
  end %% number of satellites
end %% time step

%   lonScaled(:,end)     = lonScaled(:,end-1);
%   latScaled(:,end)     = latScaled(:,end-1);
%   radScaled(:,end)     = radScaled(:,end-1);
% 
%   lon(:,end)     = lon(:,end-1);
%   lat(:,end)     = lat(:,end-1);
%   rad(:,end)     = rad(:,end-1);
  if plotLatLonIn2D
    figure
    for i=1:ns+1
      plot(lon(i,2:end-1),lat(i,2:end-1),'.');hold on;
    end
    legend
    figure
    for i=1:ns+1
      plot(lonScaled(i,2:end-1),latScaled(i,2:end-1),'.');hold on;
    end
    legend
  end
  
if writeLLRRPYData %% write files for LLR & RPY of each satellite and the reference (satellite). Latter is numbered as sat0.
  for i=1:ns+1
    rpyFile = strcat(rpyFolderPath,filesep,'sat',num2str(i-1),'_LLR_RPY.csv');
    writematrix([vizTime' latScaled(i,:)' lonScaled(i,:)' radScaled(i,:)' rollVizTime(i,:)' pitchVizTime(i,:)' yawVizTime(i,:)'],rpyFile);
  end
end

if runSPfinder %% GNSS reflectometry compuptation and visualization
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
    [timeGNSS(i,:),latGNSS(i,:),lonGNSS(i,:),radGNSS(i,:)]=kepler(cosmosTime,keplerStepSize,inclinationGNSS(i),RAANGNSS(i),v0GNSS(i),altitudeGNSS(i),radiusOfEarth);
    altGNSS(i,:)=radGNSS(i,:)-radiusOfEarth/1000;
  end
  fprintf('\nsetting up GNSS constellation time: %s seconds.',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
  
  %% compute SP location per each cubesat and each GNSS satellite
  latSP=zeros(ns+1,noOfGNSSsats,size(vizTime,2));
  lonSP=zeros(ns+1,noOfGNSSsats,size(vizTime,2));
  for i=2:ns+1
    for j=1:noOfGNSSsats
      [latSP(i,j,:), lonSP(i,j,:)]=computeSPlocation(vizTime,lat(i,:),lon(i,:),rad(i,:),timeGNSS(j,:),latGNSS(j,:),lonGNSS(j,:),altGNSS(j,:)+radiusOfEarth/1000,radiusOfEarth/1000);
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
      plotm(lat(i,:),lon(i,:),rad(i,:)-radiusOfEarth/1000,'Color',[1 0  0]);hold on;
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
    
    
%     %% plot only one line
%     GNSSsat=1;CubeSat=2;
%     plotm(latGNSS(GNSSsat,:),lonGNSS(GNSSsat,:),altGNSS(GNSSsat,:),'Color',[0 1 0]);hold on;
%     varargin={'EdgeColor',SPcolor(GNSSsat,:)};
%     circlem(squeeze(latSP(CubeSat,GNSSsat,:))',squeeze(lonSP(CubeSat,GNSSsat,:))',sizeOfSpecularPoint*ones(size(lonSP,3),1),varargin);hold on;

%% plot all lines
    for i=2:ns+1
      %% plot subsatellite points
      plotm(lat(i,:),lon(i,:),'.');hold on;

      for j=1:noOfGNSSsats
        if j<=noOfGNSSsatsArray(1)
          varargin={'EdgeColor',[1 0 0]};
        else 
          varargin={'EdgeColor',[0 0 0]};
        end
        %varargin={'EdgeColor',SPcolor(j,:)};
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
  fprintf('\n Total GNSS computation time: %s seconds.\n',num2str(posixtime(datetime('now')) - GNSScpuStartTime));
end %% end GNSS-R visualization

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

function [latSP, lonSP]=computeSPlocation(vizTime,latRec,lonRec,radRec,timeGNSS,latGNSS,lonGNSS,radGNSS,rC)
%% based on:
%% A New Approach to Determine the Specular Point of Forward Reflected GNSS Signals
%% Benjamin John Southwell Dolby Laboratories, Inc. Andrew G Dempster UNSW Sydney
%%%% input variables:
%% vizgridtime
%% latRec
%% lonRec
%% radRec
%% timeGNSS
%% latGNSS
%% lonGNSS
%% radGNSS
%% radiusOfEarth
%%%% output variables:
%% latSP
%% lonSP
  %% interpolate latGNSS/lonGNSS on lat/lon/rad grid
  latGNSSnewGrid=interp1(timeGNSS,latGNSS,vizTime,'spline',0);
  lonGNSSnewGrid=interp1(timeGNSS,lonGNSS,vizTime,'spline',0);
  radGNSSnewGrid=interp1(timeGNSS,radGNSS,vizTime,'spline',0);
  latSP=zeros(1,size(vizTime,1));
  lonSP=zeros(1,size(vizTime,1));
  
  for i=1:size(latRec,2) %% compute betas for each lat data point
    T=radGNSSnewGrid(i);
    R=radRec(i);
    %% compute thetaR, thetaT
    if latGNSSnewGrid(i)> 270 && latRec(i)<90
      THETAlat=latGNSSnewGrid(i)-latRec(i)-360;
    elseif latGNSSnewGrid(i) < 90 && latRec(i) > 270
      THETAlat=latGNSSnewGrid(i)-latRec(i)+360;
    else
       THETAlat=latGNSSnewGrid(i)-latRec(i);
    end
    
    %fprintf('\n latGNSS %f latCS %f THETAlat %f T %f R %f  ',latGNSSnewGrid(i),latRec(i),THETAlat, T, R)
    myfunLat = @(betaLat) asind(rC/T*sind(betaLat))+asind(rC/R*sind(betaLat))+THETAlat+2*betaLat-360;
    funLat = @(betaLat) myfunLat(betaLat);
    betaLat = fzero(funLat,180-THETAlat);
    betaLatOld=betaLat;
    thetaSPLat(i)=180-betaLat-asind(rC/T*sind(betaLat));
    latSPTemp(i)=latGNSSnewGrid(i)-thetaSPLat(i);    
    %fprintf('\n thetaSPLat %f betaLat %f latSPTemp %f: ',thetaSPLat(i),betaLat,latSPTemp(i))

    if lonGNSSnewGrid(i)> 270 && lonRec(i)<90
      THETAlon=lonGNSSnewGrid(i)-lonRec(i)-360;
    elseif lonGNSSnewGrid(i) < 90 && lonRec(i) > 270
      THETAlon=lonGNSSnewGrid(i)-lonRec(i)+360;
    else
      THETAlon=lonGNSSnewGrid(i)-lonRec(i);
    end
      
    %fprintf('\n lon1 %f %f %f: ',lonGNSSnewGrid(i),lonRec(i),THETAlon)
    myfunLon = @(betaLon) asind(rC/T*sind(betaLon))+asind(rC/R*sind(betaLon))+THETAlon+2*betaLon-360;
    funLon = @(betaLon) myfunLon(betaLon);
    betaLon = fzero(funLon,180-THETAlon);
    betaLonOld=betaLon;
    thetaSPLon(i)=180-betaLon-asind(rC/T*sind(betaLon));
    lonSPTemp(i)=lonGNSSnewGrid(i)-thetaSPLon(i);  
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

function [vizTime,lat,lon,rad]=kepler(cosmosTime,keplerStepSize,inclination,RAAN,v0,altitude,radiusOfEarth)

%% this function adds the global movement of the satellite based on Kepler's laws
%this function works with km instead of m %!harmonize

mu=3.986004418E14;                  %% [m3?s?2] gravitational constant
RE = radiusOfEarth/1000;          % Earth's radius                            [km]
muE = mu/1e9;    % Earth gravitational parameter             [km^3/sec^2]
wE = (2*pi/86164);  % Earth rotation velocity aorund z-axis     [rad/sec]


%% ORBIT INPUT
w=0;    %%w       = input(' Argument of perigee                  [  0,360[    w      [deg] = ');
sstTemp=RE+altitude/1000;%a       = input(' Major semi-axis                       (>6378)     a      [km]  = ');
ecc_max = sprintf('%6.4f',1-RE/sstTemp);     % maximum value of eccentricity allowed
e=0;%e       = input([' Eccentricity                         (<',ecc_max,')    e            = ']);

RAAN  = RAAN*pi/180;        % RAAN                          [rad]
w     = w*pi/180;           % Argument of perigee           [rad]
v0    = v0*pi/180;          % True anomaly at the departure [rad]
inclination     = inclination*pi/180;           % inclination                   [rad]


%% ORBIT COMPUTATION
rp = sstTemp*(1-e);               % radius of perigee             [km]
ra = sstTemp*(1+e);               % radius of apogee              [km]
Vp = sqrt(muE*(2/rp-1/sstTemp));  % velocity at the perigee       [km/s]
Va = sqrt(muE*(2/ra-1/sstTemp));  % velocity at the  apogee       [km/s]
n  = sqrt(muE./sstTemp^3);        % mean motion                   [rad/s]
p  = sstTemp*(1-e^2);             % semi-latus rectus             [km]
T  = 2*pi/n;                % period                        [s]
h  = sqrt(p*muE);           % moment of the momentum        [km^2/s]
h1 = sin(inclination)*sin(RAAN);      % x-component of unit vector h
h2 = -sin(inclination)*cos(RAAN);     % y-component of unit vector h
h3 = cos(inclination);                % z-component of unit vector h
n1 = -h2/(sqrt(h1^2+h2^2)); % x-component of nodes' line
n2 =  h1/(sqrt(h1^2+h2^2)); % y-component of nodes' line
n3 = 0;                     % z-component of nodes' line
N  = [n1,n2,n3];            % nodes' line (unit vector)
%% PRINT SOME DATAS
hours   = floor(T/3600);                   % hours   of the orbital period
minutes = floor((T-hours*3600)/60);        % minutes of the orbital period
seconds = floor(T-hours*3600-minutes*60);  % seconds of the orbital period
t0   = cosmosTime(1,1);                           % initial time          [s]
tf=cosmosTime(end,1);                             % final time
vizTime    = t0:keplerStepSize:tf;                          % vector of time        [s]
%% DETERMINATION OF THE DYNAMICS
cosE0 = (e+cos(v0))./(1+e.*cos(v0));               % cosine of initial eccentric anomaly
sinE0 = (sqrt(1-e^2).*sin(v0))./(1+e.*cos(v0));    %   sine of initial eccentric anomaly
E0 = atan2(sinE0,cosE0);                           % initial eccentric anomaly              [rad]
if (E0<0)                                          % E0€[0,2pi]
  E0=E0+2*pi;
end
tp = (-E0+e.*sin(E0))./n+t0;                       % pass time at the perigee               [s]
M  = n.*(vizTime-tp);                                    % mean anomaly                           [rad]
%% Mk = Ek - e*sin(Ek);
% Eccentric anomaly (must be solved iteratively for Ek)
E = zeros(size(vizTime,2),1);
for j=1:size(vizTime,2)
  E(j) = anom_ecc(M(j),e);                     % eccentric anomaly         [rad]
end
%% True anomaly, Argument of latitude, Radius
sin_v = (sqrt(1-e.^2).*sin(E))./(1-e.*cos(E));   % sine of true anomaly
cos_v = (cos(E)-e)./(1-e.*cos(E));               % cosine of true anomaly
v     = atan2(sin_v,cos_v);                      % true anomaly              [rad]
theta = v + w;                                   % argument of latitude      [rad]
r     = (sstTemp.*(1-e.^2))./(1+e.*cos(v));            % radius                    [km]
%% Satellite coordinates
% "Inertial" reference system ECI (Earth Centered Inertial)
xp = r.*cos(theta);                          % In-plane x position (node direction)             [km]
yp = r.*sin(theta);                          % In-plane y position (perpendicular node direct.) [km]
xs = xp.*cos(RAAN)-yp.*cos(inclination).*sin(RAAN);    % ECI x-coordinate SAT                             [km]
ys = xp.*sin(RAAN)+yp.*cos(inclination).*cos(RAAN);    % ECI y-coordinate SAT                             [km]
zs = yp.*sin(inclination);                             % ECI z-coordinate SAT                             [km]
rs = p./(1+e.*cos(theta-w));                 % norm of radius SAT                               [km]
%% GREENWICH HOUR ANGLE
%disp(' From ephemeridis you can reach the greenwich hour angle at the epoch and reset it from Aries'' point');
%disp(' Ephemeridis Almanac: http://www.nauticalalmanac.it/it/astronomia-nautica/effemeridi-nautiche.html ');
% Greenwich hour angle è l'orientamento della Terra relativo al punto di
% Ariete in un'epoca precisa. Tale valore è tabulato per ogni giorno
% dell'anno, ogni ora, minuto ecc. nelle effemeridi che si possono trovare
% negli almanacchi nautici in generale. Se non lo si conosce, o si intende
% fare un'analisi approssimativa  è consigliabile immettere 0 in modo da
% far coincidere il punto di ariete con Greenwich all'epoca.
%greenwich0 = input(' Insert Greenwich longitude respect to the vernal axis. GHA [deg] = ');
greenwich0=0;
%% SUB-SATELLITE-POINT
greenwich0 = greenwich0*pi/180;                 % Greenwich hour angle at the initial time    [rad]
rot_earth  = wE.*(vizTime-t0)+greenwich0;             % Greenwich hour angle at the time t          [rad]
for j=1:size(vizTime,2)
  if rot_earth(j) < (-pi)
    nn = ceil(rot_earth(j)/(-2*pi));
    rot_earth(j) = rot_earth(j) + nn*2*pi;
  elseif rot_earth(j) > (pi)
    nn = fix(rot_earth(j)/(2*pi));
    rot_earth(j) = rot_earth(j) - nn*2*pi;
  end
end

lat      = asin(sin(inclination).*sin(theta))/pi*180;           % Latitude             [deg]
lon      = wrapTo360((atan2(ys./rs,xs./rs)-rot_earth')/pi*180); % Longitude            [deg]
rad      = rs;                                                  % radius                [km]

end %% kepler functino

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

function [h,circlelat,circlelon] = circlem(lat,lon,radius,varargin,circleColor)
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
%%  Plot and format circle(s): 
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

function [E] = anom_ecc(M,e)
  % function [E] = anom_ecc(M,e) 
  % Risoluzione numerica dell'equazione: E-e*sin(E)=M
  % E = anomalia eccentrica [rad]
  % e = eccentricità
  % M = anomalia media [rad]
  % Si devono dare in input alla funzione due valori scalari,
  % rispettivamente M in rad ed e.
  % Il programma restituisce E [rad] con un errore di 1e-10
  % N.B. Per migliorare l'accuratezza accedere all'editor e modificare il
  % valore della variabile err, che rappresenta l'errore commesso.
  
  format long
  %x = 0;
  %sx = 1;
  %dymax = -1;
  %trovato = true;
  %while (trovato)
  %if (sx<0.2)
  % sx = 0.1 - (x/1000);
  % else
  % sx  = M-x;
  % dx  = M+x;
  % dy  = abs(cos(dx));
  %  dy2 = abs(cos(sx));
  %  end
  %  if (dymax<dy || dymax<dy2)
  %  if (dy<dy2)
  %  dymax = dy2;
  %  else
  %  dymax = dy;
  %  dy = dymax;
  %  end
  %  end 
  %  f0 = sx-e.*sin(sx)-M;
  %  f1 = dx-e.*sin(dx)-M;
  %  trovato = ((f0.*f1)>0);
  %  x = x + 0.1;
  %  end
  E = M;
  k = 1;
  err = 1e-10;
  % stabilito intervallo di ricerca (sx,dx) e max valore derivata dymax;
  while (k>err)
    y = e*sin(E)+M;
    k = abs(abs(E)-abs(y));
    E = y;
  end
  % trovato E con un errore inferiore a 0.1*10^(-4);
  %fprintf(' La soluzione E è stata trovata nell''intervallo [%2.5f,%2.5f]',sx,dx);
  %fprintf(' \n errore inferiore a %1.2e: [rad] E = %3.5f',err,E);
end
