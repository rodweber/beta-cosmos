function ECEFprocessing(this,vizScale,ns,altitude,radiusOfEarth)
%% input:
% ns: number of satellites
% altitude: altitude used for visualization. The actual altitude may change for long periods of
% operations. However, this will be so small that it is irrelevant for visualization
%% output:
% nonef
% however, files are written that are used by cosmosVIZ

%% switches
plotLatLonIn2D=0;
writeLLRRPYData=1;

%% global parameter for formation
%inclination=0; %% [-90, 90] [deg]  %% equatorial
%inclination=35; %% [-90, 90] [deg] %% CYGNSS
%inclination=52; %% [-90, 90] [deg] %% ISS
%inclination=90; %% [-90, 90] [deg] %% polar
inclination=97; %% [-90, 90] [deg] %% (approximately) SSO

RAAN=0; %%RAAN    = input(' Right Ascension of Ascendent Node    [  0,360[    RAAN   [deg] = ');
%v0=225;   %%v0      = input(' True anomaly at the departure        [  0,360[    v0     [deg] = ');
v0=0;   %%v0      = input(' True anomaly at the departure        [  0,360[    v0     [deg] = ');
%% this is the timestep used in the following, also for the GNSSR processing.
keplerStepSize=1; %% [s]


fprintf('\nECEF processing...');


%% read data from telemery files
% JT: this works only if the telemetry data of all satellites is equal in size 
% and sync'ed. Sooner or later, this needs to become more versatile allowing individual times for each satellite
for i=1:ns
  tempTime=readmatrix(strcat('TimeVectorTM',num2str(i),'.csv'));  
  tempSatStates=readmatrix(strcat('SatStatesTM',num2str(i),'.csv')); 
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
[vizTime,latitude,longitude,radius] = this.keplerPropagation(cosmosTime,keplerStepSize,inclination,RAAN,v0,altitude,radiusOfEarth);

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

%% plot latitude and longitude in X-Y plot
if plotLatLonIn2D
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
  for i=1:ns+1
    writematrix([vizTime' latScaled(i,:)' lonScaled(i,:)' radScaled(i,:)' rollVizTime(i,:)' pitchVizTime(i,:)' yawVizTime(i,:)'  ],strcat('sat',num2str(i-1),'_LLR_RPY_Scaled.csv'));
    writematrix([vizTime' lat(i,:)' lon(i,:)' rad(i,:)' ],strcat('sat',num2str(i-1),'_LLR.csv'));
  end
  fprintf('\nECEF files...written');
end

fprintf('\nECEF processing...done');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
