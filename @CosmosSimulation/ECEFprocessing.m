function ECEFprocessing(this, altitude, inclination, orbitLAN, vizScale, ...
    keplerStepSize, v0, plotLatLonIn2D, writeLLRRPYData, parentCoordFolder,...
    llrNotScaledFolderName, llrScaledFolderName, xyzScaledFolderName, enablePlot)
%% input:
% altitude: altitude used for visualization. The actual altitude may change for long periods of
%   operations. However, this will be so small that it is irrelevant for visualization
% inclination, [-90, 90] [deg] %% equatorial, 35 CYGNSS, 52 ISS, 90 polar, 97 (approximately) SSO
% keplerStepSize, [s] %% timestep used in the following, also for the GNSSR processing
% RAAN, [0,360[ [deg] %% Right Ascension of Ascendent Node
% v0,[  0,360[,[deg] %% True anomaly at the departure
% plotLatLonIn2D %% switch
% writeLLRRPYData %% switch
%% output:
% nonef
% however, files are written that are used by cosmosVIZ

%% Start
fprintf('\nECEF processing ... ');

% Set some parameters.
ns = this.NumSatellites;
radiusOfEarth = this.Orbits(1).MeanEarthRadius;

% Get path to telemetry files from CosmosSimulation object.
tmFolderPath = this.TelemetryPath;

% Set path for coordinate folders.
llrNotScaledFolderPath = strcat(parentCoordFolder,filesep,llrNotScaledFolderName);
llrScaledFolderPath = strcat(parentCoordFolder,filesep,llrScaledFolderName);
xyzScaledFolderPath = strcat(parentCoordFolder,filesep,xyzScaledFolderName);

% Create required folders and delete old files.
[~,~,~] = mkdir(parentCoordFolder); % [status, msg, msgID]
[~,~,~] = mkdir(llrNotScaledFolderPath); % [status, msg, msgID]
[~,~,~] = mkdir(llrScaledFolderPath); % [status, msg, msgID]
[~,~,~] = mkdir(xyzScaledFolderPath); % [status, msg, msgID]
delete(strcat(parentCoordFolder,filesep,'*.csv'));
delete(strcat(llrNotScaledFolderPath,filesep,'*.csv'));
delete(strcat(llrScaledFolderPath,filesep,'*.csv'));
delete(strcat(xyzScaledFolderPath,filesep,'*.csv'));

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
      sstX =zeros(ns,timeSteps);
      sstY =zeros(ns,timeSteps);
      sstZ =zeros(ns,timeSteps);
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

%% Enable Rod method for circular orbits and all interdependent sections
ENABLE_ROD_METHOD = true;
% Vector of time.
t0 = cosmosTime(1,1); % Initial time [s].
tf = cosmosTime(end,1); % Final time [s].
vizTime = (t0:keplerStepSize:tf)'; % [s] [TIMEx1 vector]
dataLength = length(vizTime);
% Relative positions at the timestamps of cosmos simulation.
temp_xlocal = sstX'; % [TIMExNUMSATS vector]
temp_ylocal = sstY'; % [TIMExNUMSATS vector]
temp_zlocal = sstZ'; % [TIMExNUMSATS vector]
% Interpolate relative positions at the visualization timestamps.
xlocal = zeros(dataLength,ns);
ylocal = zeros(dataLength,ns);
zlocal = zeros(dataLength,ns);
for n = 1:ns
  xlocal(:,n) = interp1(cosmosTime(:,n), temp_xlocal(:,n), vizTime);
  ylocal(:,n) = interp1(cosmosTime(:,n), temp_ylocal(:,n), vizTime);
  zlocal(:,n) = interp1(cosmosTime(:,n), temp_zlocal(:,n), vizTime);
end

%% rotate relative coordinate system to align x+ with flight direction (inclination)
for i=1:ns
  for j=1:size(sstX,2)
    sstTemp=rotz(inclination-90)*[sstX(i,j) sstY(i,j) 0]';
  end
  sstX(i,j)=sstTemp(1);
  sstY(i,j)=sstTemp(2);
end

%% compute orbit
[vizTime,latitude,longitude,radius] = this.keplerPropagation(cosmosTime,keplerStepSize,inclination,orbitLAN,v0,altitude,radiusOfEarth);
vizTime = vizTime'; % transform to TIMEx1 column vector
if(ENABLE_ROD_METHOD)
  [vizTime,x,y,z] = this.getPosECEF(cosmosTime,keplerStepSize,inclination,orbitLAN,v0,altitude,radiusOfEarth,...
    xlocal,ylocal,zlocal,vizScale,ns);
end

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

%{
%% off set of the formation satellites
for i=1:length(vizTime)
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
for i=1:length(vizTime)
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
%}
fprintf('done.\n');

%% plot latitude and longitude in X-Y plot
if plotLatLonIn2D && enablePlot
    figure
    for i=1:ns+1
      plot(lon(i,2:end-1),lat(i,2:end-1),'.');hold on;
    end
    title('latlon');
    legend
    figure
    for i=1:ns+1
      plot(lonScaled(i,2:end-1),latScaled(i,2:end-1),'.');hold on;
    end
    title('latlonscaled');
    legend
end

%% write files for LLR & RPY of each satellite and the reference (satellite). Latter is numbered as sat0.
if writeLLRRPYData 
  % Create inclination vector to write to output file.
  incVec = inclination * ones(size(vizTime));
  lanVec = orbitLAN * ones(size(vizTime));
  fprintf('Writing ECEF files ... ');
  for i=1:ns+1
    llrFileName = strcat(llrNotScaledFolderPath,filesep,'sat',num2str(i-1),'_LLR.csv');
    llrrpyScaledFileName = strcat(llrScaledFolderPath,filesep,'sat',num2str(i-1),'_LLR_RPY_Scaled.csv');
    xyzrpyScaledFileName = strcat(xyzScaledFolderPath,filesep,'sat',num2str(i-1),'_XYZ_RPY_Scaled.csv');
    if(ENABLE_ROD_METHOD)
      writematrix([vizTime x(:,i) y(:,i) z(:,i) rollVizTime(i,:)' pitchVizTime(i,:)' yawVizTime(i,:)' incVec lanVec], xyzrpyScaledFileName);
    end
    %writematrix([vizTime latScaled(i,:)' lonScaled(i,:)' radScaled(i,:)' rollVizTime(i,:)' pitchVizTime(i,:)' yawVizTime(i,:)' incVec], llrrpyScaledFileName);
    %writematrix([vizTime lat(i,:)' lon(i,:)' rad(i,:)'], llrFileName);
  end
  fprintf('done.\n');
end

end
