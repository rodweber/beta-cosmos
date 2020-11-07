function visualizationLONLATALT(vizScale,ns,altitude)
%% input:
% ns: number of satellites
% altitude: altitude used for visualization. The actual altitude may change for long periods of
% operations. However, this will be so small that it is irrelevant for visualization
%% output:
% none
% however, files are written that are used by cosmosVIZ

%% global parameter for formation
inclination=90; %% [-90, 90] [deg]
radiusOfEarth=6371000;          %% [m]
RAAN=0; %%RAAN    = input(' Right Ascension of Ascendent Node    [  0,360[    RAAN   [deg] = ');
v0=0;   %%v0      = input(' True anomaly at the departure        [  0,360[    v0     [deg] = ');


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
  sstX(i,:)=tempSatStates(:,1)';
  sstY(i,:)=tempSatStates(:,2)';
  sstZ(i,:)=tempSatStates(:,3)';
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
  [vizTime,latitude,longitude,radius]=kepler(cosmosTime,inclination,RAAN,v0,altitude,radiusOfEarth);
  
  %% interpolate relative position on visualization time steps
  for i=1:ns
    %% interpolate on viz-time-grid
    sstXvizTime(i,:)=interp1(cosmosTime(:,i),sstX(i,:),vizTime);
    sstYvizTime(i,:)=interp1(cosmosTime(:,i),sstY(i,:),vizTime);
    sstZvizTime(i,:)=interp1(cosmosTime(:,i),sstZ(i,:),vizTime);    

    %% interpolate on viz-time-grid and apply visualization scaling
    sstXvizTimeScaled(i,:)=interp1(cosmosTime(:,i),sstX(i,:)*vizScale,vizTime);
    sstYvizTimescaled(i,:)=interp1(cosmosTime(:,i),sstY(i,:)*vizScale,vizTime);
    sstZvizTimeScaled(i,:)=interp1(cosmosTime(:,i),sstZ(i,:)*vizScale,vizTime);
    
    %% interpolate Euler angles
    rollVizTime(i,:)  =interp1(cosmosTime(:,i),squeeze(roll(i,:)),vizTime);
    pitchVizTime(i,:) =interp1(cosmosTime(:,i),squeeze(pitch(i,:)),vizTime);
    yawVizTime(i,:)   =interp1(cosmosTime(:,i),squeeze(yaw(i,:)),vizTime);
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
  rollVizTime   = [zeros(1,size(pitchVizTime,2)); rollVizTime];
  pitchVizTime  = [zeros(1,size(pitchVizTime,2)); pitchVizTime];
  yawVizTime    = [zeros(1,size(pitchVizTime,2)); yawVizTime];
  
  %% off set of the formation satellites
  for i=1:size(vizTime,2)-1
    for j=1:ns
      localTransformedCoSystemScaled = rotz(inclination*180/pi-90)*[sstXvizTimeScaled(j,i) ; sstYvizTimescaled(j,i); sstZvizTimeScaled(j,i)];
      offPlaneOffsetAngleScaled(j,i) = asind(  localTransformedCoSystemScaled(1)/1000 / (radScaled(1,i) +  localTransformedCoSystemScaled(3)/1000) );
      inPlaneOffsetAngleScaled(j,i)  = asind(  localTransformedCoSystemScaled(2)/1000 / (radScaled(1,i)  + localTransformedCoSystemScaled(3)/1000) );

      localTransformedCoSystem = rotz(inclination*180/pi-90)*[sstXvizTime(j,i) ; sstYvizTime(j,i); sstZvizTime(j,i)];
      offPlaneOffsetAngle(j,i) = asind(  localTransformedCoSystem(1)/1000 / (rad(1,i) +  localTransformedCoSystem(3)/1000) );
      inPlaneOffsetAngle(j,i)  = asind(  localTransformedCoSystem(2)/1000 / (rad(1,i)  + localTransformedCoSystem(3)/1000) );
    end
  end
  
  %% Lat, Lon, Rad of formation satellites 
  for i=1:size(vizTime,2)-1
    for j=2:ns+1

      lonScaled(j,i)     = wrapTo360(lonScaled(1,i)+offPlaneOffsetAngleScaled(j-1,i));         % Longitude            [deg]
      latScaled(j,i)     = latScaled(1,i)+inPlaneOffsetAngleScaled(j-1,i);                     % Latitude             [deg]
      radScaled(j,i)     = radScaled(1,i)'+sstZvizTimeScaled(j-1,i)/1000;                      % radius                [km]

      lon(j,i)     = wrapTo360(lon(1,i)+offPlaneOffsetAngle(j-1,i));         % Longitude            [deg]
      lat(j,i)     = lat(1,i)+inPlaneOffsetAngle(j-1,i);                     % Latitude             [deg]
      rad(j,i)     = rad(1,i)'+sstZvizTime(j-1,i)/1000;                      % radius                [km]

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

  lonScaled(:,end)     = lonScaled(:,end-1);
  latScaled(:,end)     = latScaled(:,end-1);
  radScaled(:,end)     = radScaled(:,end-1);

  lon(:,end)     = lon(:,end-1);
  lat(:,end)     = lat(:,end-1);
  rad(:,end)     = rad(:,end-1);
       
  %% write files for LLR & RPY of each satellite and the reference (satellite). Latter is numbered as sat0.
  for i=1:ns+1
    writematrix([vizTime' latScaled(i,:)' lonScaled(i,:)' radScaled(i,:)' rollVizTime(i,:)' pitchVizTime(i,:)' yawVizTime(i,:)'  ],strcat('sat',num2str(i-1),'_LLR_RPY.csv'));
  end
  
  if 1 %% GNSS reflectometry visualization
    %% compute or define GNSS satellites position
    noOfGNSSsats=2;
    
    inclinationGNSS(1)=60;
    RAANGNSS(1)=0;
    v0GNSS(1)=0;
    altitudeGNSS(1)=20e6;
    
    inclinationGNSS(2)=30;
    RAANGNSS(2)=0;
    v0GNSS(2)=0;
    altitudeGNSS(2)=20e6;

    for i=1:noOfGNSSsats
     [timeGNSS(i,:),latGNSS(i,:),lonGNSS(i,:),radGNSS(i,:)]=kepler(cosmosTime,inclinationGNSS(i),RAANGNSS(i),v0GNSS(i),altitudeGNSS(i),radiusOfEarth);
     altGNSS(i,:)=radGNSS(i,:)-radiusOfEarth/1000;
    end
    %% compute SP location per each cubesat and each GNSS satellite
    for i=2:ns+1
      for j=1:noOfGNSSsats
        [latSP(i,j,:), lonSP(i,j,:)]=computeSPlocation(vizTime,lat(i,:),lon(i,:),rad(i,:),timeGNSS(j,:),latGNSS(j,:),lonGNSS(j,:),altGNSS(j,:)+radiusOfEarth/1000,radiusOfEarth/1000);
      end
    end
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
          %plotm(squeeze(latSP(i,j,:))',squeeze(lonSP(i,j,:))',ones(1,size(lonSP,3)),'o','Color',[0 0 0]);hold on;
          plotm(squeeze(latSP(i,j,:))',squeeze(lonSP(i,j,:))',ones(1,size(lonSP,3)),'o');hold on;
        end
      end
      view(90,0)
       
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
        for j=1:noOfGNSSsats
          plotm(squeeze(latSP(i,j,:))',squeeze(lonSP(i,j,:))','o');hold on;
        end
      end       
       
  end %% end GNSS-R visualization
  
end




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
      latSP(i)=latSPTemp(i);
      lonSP(i)=lonSPTemp(i);
    else
      latSP(i)=0;
      lonSP(i)=0;      
    end
  %  if lonSP(i)>90
  %    input('x');
  %  end
  end
 % latSP
 % lonSP
 % input('xXXX');
  
end


function [vizTime,lat,lon,rad]=kepler(cosmosTime,inclination,RAAN,v0,altitude,radiusOfEarth)

%% this function adds the global movement of the satellite based on Kepler's laws
%this function works with km instead of m %!harmonize

mu=3.986004418E14;                  %% [m3?s?2] gravitational constant
%dimsst = size(sstx);
%sstx = reshape(sstx,[dimsst(2:end) 1]);
%ssty = reshape(ssty,[dimsst(2:end) 1]);
%sstz = reshape(sstz,[dimsst(2:end) 1]);

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
vis_step=2;       %step = input(' Time step.        [s] step = ');  % time step             [s]
vizTime    = t0:vis_step:tf;                          % vector of time        [s]
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


end
