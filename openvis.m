%% File to open visualization for Cosmos Beta in MATLAB Simulink

%% Constants
% Earth rotation velocity around Z-axis.
EARTH_ROT = (2*pi/86164); % [rad/s]

%% Set paths
% The parameters below normally should never change.
warning on verbose;

% Inform the name of this file without the extension "m".
THIS_FILE_NAME = 'openvis';

% Inform the name of the project file for the Cosmos simulation.
PROJECT_FILE_NAME = 'cosmosVisualization'; % Without extension 'prj'.

if(~isdeployed)
  % Get directory path of the active file in MATLAB's Editor.
  [pathActiveFile,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
  addpath(genpath(pathActiveFile)); % Append path to the current MATLAB path.
  
  % Get path for directory of the file name set in THIS_FILE_NAME.
  [pathCosmos,~,~] = fileparts(which(THIS_FILE_NAME));
  addpath(genpath(pathCosmos)); % Append path to the current MATLAB path.
  
  % Add path to ancillary folders.
  addpath(strcat(pathCosmos,filesep,'config'));
  addpath(strcat(pathCosmos,filesep,'utils'));
  
  % Get path for directory of the simulation project.
  [pathVisualization,~,~] = fileparts(which([PROJECT_FILE_NAME,'.prj']));
  
  % Change working directory to the directory of this m-file.
  cd(pathCosmos);
end

% Check if there is any Matlab project already open.
if(isempty(matlab.project.rootProject))
  % If there is no project open, launch project from path.
  proj = openProject(pathVisualization);
  fprintf('\nTo prevent issues, check if project path is correct:\n');
  fprintf('%s\n\n',proj.RootFolder);
else
  % If there is a project already open, get the Project object.
  proj = currentProject();
  % Check if the name of the project is correct.
  if(~strcmp(proj.Name,PROJECT_FILE_NAME))
    % If name is different, show error.
    error(['\nProject ''%s'' is already open. ',...
      'Close it before running ''%s''.'], ...
      proj.Name, PROJECT_FILE_NAME);
  else
    % If name is correct, show message.
    fprintf(2,'\nAttention: Project with name ''%s'' already open.\n',proj.Name);
    fprintf('To prevent issues, check if project path is correct:\n');
    fprintf('%s\n\n',proj.RootFolder);
  end
end

% Set timer start time.
timeScriptStart = posixtime(datetime('now')); % Posixtime [seconds].


%% Read parameters from JSON files
vis = readjson('configVisualization.json');

% Set parameter to automatically run Simulink visualization.
AUTORUN = vis.AutoRun; % [true: 1 | false: 0]

% Set path for coordinates folder.
% IMPORTANT: The FIRST file in the coordinates folder MUST BE the reference for
% the other satellites.
COORD_FOLDER = strcat(vis.ParentCoordFolder,filesep,vis.CoordFolder);
COORD_FOLDER = fullfile(pathCosmos,COORD_FOLDER);

% Set struct for selected camera mode.
CAM_MODE = struct('time',0,'signals',struct('dimensions',0,'values',0));

% Set parameter to automatically smooth changes in satellite orientations.
SMOOTH_ENABLE = vis.SmoothSatOrientationChanges; % [true: 1 | false: 0]

% If smoothing is enabled, set num of data points for computing smoothed value.
SMOOTH_SPAN = vis.SmoothingDataSpan;
SMOOTH_METHOD_LIST = vis.SmoothingMethodList;
SMOOTH_METHOD = SMOOTH_METHOD_LIST{vis.SmoothingMethodChosen};

% Option to override lat-long-radius values for testing purposes.
OVERRIDE_LLR = vis.OverrideLLR;

% Option to override roll-pitch-yaw values for testing purposes.
OVERRIDE_RPY = vis.OverrideRPY;


%% Prepare coordinate files
disp('Satellite coordinate files:');
disp(strcat(pathCosmos,filesep));
fprintf('  %s\n',strcat(vis.ParentCoordFolder,filesep,vis.CoordFolder));
% Read names of coordinate files.
coordfiles = {dir(strcat(COORD_FOLDER,filesep,'*.csv')).name};
% Rearrange files so that the local coordinate is the last one in the list.
% This step is necessary for correct visualization in the 3D simulink world.
coordfiles = circshift(coordfiles,-1);
disp(coordfiles');

% Clear satellite struct.
clear simsat;

% Define satellite struct.
simsat(1).time  = -1;

simsat(1).x      = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).y      = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).z      = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).xyz    = struct('time',0,'signals',struct('dimensions',0,'values',0));

simsat(1).lat    = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).long   = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).sma    = struct('time',0,'signals',struct('dimensions',0,'values',0));

simsat(1).pos    = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).posU1  = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).posU2  = struct('time',0,'signals',struct('dimensions',0,'values',0));

simsat(1).pitch  = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).yaw    = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).roll   = struct('time',0,'signals',struct('dimensions',0,'values',0));

simsat(1).rot    = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).inc    = -1;
simsat(1).lanDeg = -1;
simsat(1).lanRad = -1;

% Allocate memory for array of satellites.
numsats = length(coordfiles);
simsat = repmat(simsat(1),numsats,1);


%% Set parallel pool
createparpool(numsats);

% Create data queue for parallel pool.
dq = parallel.pool.DataQueue;

% Define function to call when new data is received on the DataQueue.
afterEach(dq, @disp);


%% Read data from coordinate files
fprintf('%s','Reading data from coordinate files ... ');
% Set timer start time.
timeReadStart = posixtime(datetime('now')); % Posixtime [seconds].
for n = 1:numsats
  coord = readmatrix(strcat(COORD_FOLDER,filesep,coordfiles{n}));
  timestamps = coord(:,1); % [seconds]
  x = coord(:,2); % [m]
  y = coord(:,3); % [m]
  z = coord(:,4); % [m]
  pitch = coord(:,6); % [degrees]
  yaw   = coord(:,7); % [degrees]
  roll  = coord(:,5); % [degrees]
  inclinationDeg = coord(1,8); % [degrees]
  orbitLANDeg = coord(1,9); % [degrees]
  
  % Get length of data samples.
  dataLength = length(timestamps);
  % Set column array of zeros with length of data samples.
  zeroArray = zeros(dataLength,1); % [TIMEx1 vector]
  % Set column array of ones with length of data samples.
  onesArray = ones(dataLength,1); % [TIMEx1 vector]
  
  % Update camera mode struct.
  CAM_MODE.time = timestamps;
  CAM_MODE.signals.dimensions = 1;
  CAM_MODE.signals.values = vis.SelectedCameraMode*onesArray; % [TIMEx1 vector]
  
  if(OVERRIDE_LLR)
    % Values to override for each of the satellites.
    smaKm = onesArray * 6871;
    if n == 1
      latDeg  = onesArray * -20;
      longDeg = onesArray * 0;
    elseif n == 2
      latDeg  = onesArray * 0;
      longDeg = onesArray * 0;
    elseif n == 3
      latDeg  = onesArray * 20;
      longDeg = onesArray * 0;
    else
      % Reference (n == 4).
      latDeg  = onesArray * 0;
      longDeg = onesArray * 0;
    end
  end
  
  if(OVERRIDE_RPY)
    % Values to override for each of the satellites.
    if n == 1
      pitch = onesArray * 0;
      yaw   = onesArray * 0;
      roll  = onesArray * 45;
    elseif n == 2
      pitch = onesArray * 0;
      yaw   = onesArray * 90;
      roll  = onesArray * 45;
    elseif n == 3
      pitch = onesArray * 45;
      yaw   = onesArray * 0;
      roll  = onesArray * 90;
    else
      % Reference (n == 4).
    end
  end
  
  if(SMOOTH_ENABLE)
    pitch = smooth(pitch,SMOOTH_SPAN,SMOOTH_METHOD);
    yaw   = smooth(yaw  ,SMOOTH_SPAN,SMOOTH_METHOD);
    roll  = smooth(roll ,SMOOTH_SPAN,SMOOTH_METHOD);
  end
  
  % Convert data to SI units.
  pitch = pitch * (pi/180); % [rad]
  yaw   = yaw   * (pi/180); % [rad]
  roll  = roll  * (pi/180); % [rad]
  inclinationRad = inclinationDeg * (pi/180); % [rad]
  orbitLANRad = orbitLANDeg * (pi/180); % [rad]
  
  % Place data into satellite struct.
  simsat(n).x.time = timestamps;
  simsat(n).x.signals.dimensions = 1;
  simsat(n).x.signals.values = x; % [m] [TIMEx1 vector]
  
  simsat(n).y.time = timestamps;
  simsat(n).y.signals.dimensions = 1;
  simsat(n).y.signals.values = y; % [m] [TIMEx1 vector]
  
  simsat(n).z.time = timestamps;
  simsat(n).z.signals.dimensions = 1;
  simsat(n).z.signals.values = z; % [m] [TIMEx1 vector]
  
  simsat(n).xyz.time = timestamps;
  simsat(n).xyz.signals.dimensions = 3;
  simsat(n).xyz.signals.values = [x y z]; % [TIMEx3 vector]
  
  simsat(n).pos.time = timestamps;
  simsat(n).pos.signals.dimensions = 3;
  simsat(n).pos.signals.values = [zeroArray zeroArray zeroArray];
  
  simsat(n).posU1.time = timestamps;
  simsat(n).posU1.signals.dimensions = 3;
  simsat(n).posU1.signals.values = [zeroArray zeroArray zeroArray];
  
  simsat(n).posU2.time = timestamps;
  simsat(n).posU2.signals.dimensions = 3;
  simsat(n).posU2.signals.values = [zeroArray zeroArray zeroArray];
  
  simsat(n).rot.time = timestamps;
  simsat(n).rot.signals.dimensions = 4;
  simsat(n).rot.signals.values = zeros(dataLength,4); % [TIMEx4 vector]
  
  simsat(n).pitch.time = timestamps;
  simsat(n).pitch.signals.dimensions = 1;
  simsat(n).pitch.signals.values = pitch; % [rad]
  
  simsat(n).yaw.time = timestamps;
  simsat(n).yaw.signals.dimensions = 1;
  simsat(n).yaw.signals.values = yaw; % [rad]
  
  simsat(n).roll.time = timestamps;
  simsat(n).roll.signals.dimensions = 1;
  simsat(n).roll.signals.values = roll; % [rad]
  
  simsat(n).inc = inclinationRad;
  simsat(n).lanDeg = orbitLANDeg;
  simsat(n).lanRad = orbitLANRad;
  simsat(n).time = timestamps; % [seconds]
end
fprintf('%s','done. ');
% Set timer stop time.
timeReadStop = posixtime(datetime('now')); % Posixtime [seconds].
timeReadDuration = timeReadStop - timeReadStart;
fprintf('(%.2f seconds)\n',timeReadDuration);


%% Make transformations for 3D visualization
% Get last timestamp of the time vector, set it as simulation time.
stopTime = simsat(1).time(end);

% Get data length.
dataLength = length(simsat(1).time);
fprintf('%s%d.\n','Data length of each coordinate file: ',dataLength);

parsat = simsat;

% Get each time-coordinate-orientation triple for each satellite.
fprintf('%s\n','Processing data from coordinate files ... ');
spmd(numsats)
  
  n = labindex;
  
  % Set timer start time.
  timeProcStart = posixtime(datetime('now')); % Posixtime [seconds].
  msg = sprintf('%3.d/%d ... processing ...',n,numsats);
  send(dq, msg);
  for i = 1:dataLength
    
    % Get time, in seconds, since beginning of simulation.
    time = parsat(n).time(i);
    
    % Get coordinate values for the current time.
    x = parsat(n).x.signals.values(i);
    y = parsat(n).y.signals.values(i);
    z = parsat(n).z.signals.values(i);
    
    parsat(n).pos.signals.values(i,1) = x;
    parsat(n).pos.signals.values(i,2) = y;
    parsat(n).pos.signals.values(i,3) = z;
    
    % Get roll, pitch, yaw angles.
    rollAngle  = parsat(n).roll.signals.values(i);  % [rad]
    pitchAngle = parsat(n).pitch.signals.values(i); % [rad]
    yawAngle   = parsat(n).yaw.signals.values(i);   % [rad]
    
    % Compute unit vector from x, y, z.
    xyzUnit = [x;y;z] / norm([x;y;z]); % [3x1 unit vector]
    
    
    
% posVector = [x y z];
% unitPositionVector1 = posVector / norm(posVector);
%     
% % Set U1.
% parsat(n).posU1.signals.values(i,1) = unitPositionVector1(1);
% parsat(n).posU1.signals.values(i,2) = unitPositionVector1(2);
% parsat(n).posU1.signals.values(i,3) = unitPositionVector1(3);
% 
% % For each pair of lat-long, compute unit vector to center of Earth.
% % Initial vector is [SMA, 0, 0].
% % Normalize vector to start with vector0 = [1 0 0].
% unitPositionVector2 = [1 0 0];
% 
% % Compute latitude and longitude rotations.
% a = -lat;
% rotLat = [cos(a/2), 0, sin(a/2), 0]; % [quaternion]
% a = long;
% rotLong = [cos(a/2), 0, 0, sin(a/2)]; % [quaternion]
% 
% % Update unit vector pointing to satellite position.
% unitPositionVector2 = quatrotate(rotLat, unitPositionVector2);
% unitPositionVector2 = quatrotate(rotLong, unitPositionVector2);
% 
% % Set U2.
% parsat(n).posU2.signals.values(i,1) = unitPositionVector2(1);
% parsat(n).posU2.signals.values(i,2) = unitPositionVector2(2);
% parsat(n).posU2.signals.values(i,3) = unitPositionVector2(3);
    
    
    
    
    % ------------------------------------------------------------------------ %
    % Method 1 (more complex, stop here, implement Method 2)
    % For each pair of lat-long, compute vector to center of Earth.
    % Start with vector1 = [1 0 0].
    
    % Rotate on ECEF's Z-axis with angle = longitude.
    % Outputs rotated vector2.
    
    % Compute ortogonal vector between vector2 and Z-axis positive.
    % Ortogonal vector must be the rotation axis for vector2 to meet Z+.
    % Check it.
    
    % Compute rotation matrix on ortogonal vector with angle = latitude.
    % Use Rodrigues' rotation equations.
    
    % Use rotation matrix on vector2.
    % Outputs vector3.
    % ------------------------------------------------------------------------ %
    
    % Method 2: Pseudo-code here, implementation has been changed below
    % For each pair of lat-long, compute vector to center of Earth.
    % Initial vector is [SMA, 0, 0].
    % Normalize vector to start with vector0 = [1 0 0].
    
    % Rotate to set initial attitude:
    %   - payload to Nadir, and
    %   - wings perpendicular to Equator line.
    % Output vector1.
    
    % Rotate on ECEF's Y-axis with angle -latitude.
    % Output vector2.
    
    % Rotate on ECEF's Z-axis with angle +longitude.
    % Output vector3.
    
    % Compute Rodrigues' matrix for rotation around vector3.
    % Use Rodrigues' matrix with two different angles.
    % Case 1: angle = +orbit_inclination
    % Case 2: angle = -orbit_inclination
    
    % If the lat-long point is in the ascending portion of the orbit, the 
    % inclination will be positive.
    % If the lat-long point is in the descending portion of the orbit, the 
    % inclination will be negative.
    
    % Compute the derivative of the latitude.
    % If derivative is positive, orbit is ascending.
    % If derivative is negative, orbit is descending.
    
    % Therefore,
    % Derivative+ : ascending : angle = +orbit_inclination
    % Derivative- : ascending : angle = -orbit_inclination
    % ------------------------------------------------------------------------ %
    
    
    %% Correction of initial orientation
    % Rotate to set initial attitude:
    %   - payload pointing to Nadir, and
    %   - wings perpendicular to Equator line.
    % To correct orientation of the 3D satellite to initial attitude,
    % Rotate -90deg (-pi/2) on ECEF's Y-axis.
    a = -pi/2;
    rot = [cos(a/2), 0, sin(a/2), 0]; % [quaternion]
    
    
    %% Initial body axes
    % Define local satellite axes according to ECEF's XYZ axes.
    localYaw   = [1 0 0]'; % [3x1 U-vector]
    localRoll  = [0 1 0]'; % [3x1 U-vector]
    localPitch = [0 0 1]'; % [3x1 U-vector]
    
    
    %% Satellite rotation for orbital inclination
    % Positive rotation on ECEF's X-axis.
    a = parsat(n).inc; % [rad]
    rotInc = [cos(a/2), sin(a/2), 0, 0]; % [quaternion]
    rot = quatmultiply(rotInc, rot);
    % Correct local body axes.
    adeg = a * 180/pi; % [deg]
    rotInc = rotx(adeg);
    localRoll  = rotInc*localRoll;  % [3x1 U-vector]
    localPitch = rotInc*localPitch; % [3x1 U-vector]
    localYaw   = rotInc*localYaw;   % [3x1 U-vector]
    
    
    %% Total Earth rotation since beginning of simulation
    aEarthRad = EARTH_ROT * time; % [rad]
    aEarthDeg = aEarthRad * 180/pi; % [deg]
    % Correct for orbit longitude of ascending node (LAN).
    aEarthRad = aEarthRad - parsat(n).lanRad; % [rad]
    aEarthDeg = aEarthDeg - parsat(n).lanDeg; % [deg]
    
    
    %% Satellite correction for Earth rotation
    % Negative rotation on ECEF's Z-axis.
    a = -aEarthRad; % [rad]
    rotEarthFix = [cos(a/2), 0, 0, sin(a/2)]; % [quaternion]
    rot = quatmultiply(rotEarthFix, rot);
    % Correct local body axes.
    localRoll  = rotz(-aEarthDeg)*localRoll;  % [3x1 U-vector]
    localPitch = rotz(-aEarthDeg)*localPitch; % [3x1 U-vector]
    localYaw   = rotz(-aEarthDeg)*localYaw;   % [3x1 U-vector]
    
    
    %% Normal vector of the orbital plane
    % Set unit vector pointing towards ECEF's +Z-axis.
    %   unitVector = [0; 0; 1];
    % Rotate unit vector on ECEF's X-axis, with angle equal to the orbital
    % inclination. Rotation angle in degrees:
    angX = parsat(n).inc * 180/pi; % [deg]
    % Apply rotation on unit vector. The resulting vector represents the normal
    % vector of the orbital plane:
    normalOrbitVector = rotx(angX)*[0;0;1]; % [3x1 U-vector]
    % Consider retrograde motion of the orbit.
    % Rotate normal orbit vector around Z-axis.
    normalOrbitVector = rotz(-aEarthDeg)*normalOrbitVector; % [3x1 U-vector]
    
    
    %% Angle from Equator line to current orbital position
    % Equator point considering orbit retrograde motion.
    v1 = rotz(-aEarthDeg)*[1;0;0]; % [3x1 U-vector]
    % Current orbit position.
    v2 = xyzUnit; % [3x1 U-vector]
    % Normal vector of plane between Equator point and current position.
    normalVector = normalOrbitVector; % [3x1 U-vector]
    % From function a = vecangle360(v1,v2,n).
    x = cross(v1,v2);
    c = sign(dot(x,normalVector)) * norm(x);
    angPos = atan2(c,dot(v1,v2)); % [rad]
    
    
    %% Pitch correction for orbital position
    % Use corrected pitch axis as axis for pitch rotation.
    u = localPitch;
    % Rotate with an angle equal to the angle of the current orbital position.
    a = angPos; % [rad]
    % Compute Rodrigues rotation matrix.
    W = [  0  -u(3)  u(2);
         u(3)   0   -u(1);
        -u(2)  u(1)   0 ];
    I = eye(3); % identity matrix 3x3
    rodriguesRotMatrix = I + sin(a)*W + (2*sin(a/2)^2)*(W^2);
    % Convert Rodrigues rotation matrix to quaternions.
    pitchCorrection = rotm2quat(rodriguesRotMatrix);
    rot = quatmultiply(pitchCorrection, rot);
    % Correct local body axes.
    localRoll  = rodriguesRotMatrix*localRoll;  % [3x1 U-vector]
    localPitch = rodriguesRotMatrix*localPitch; % [3x1 U-vector]
    localYaw   = rodriguesRotMatrix*localYaw;   % [3x1 U-vector]
    
    
    %% Roll
    u = localRoll; % Set axis for rotation.
    a = rollAngle; % [rad]
    W = [  0  -u(3)  u(2);
         u(3)   0   -u(1);
        -u(2)  u(1)   0 ];
    I = eye(3); % identity matrix 3x3
    rodriguesRotMatrix = I + sin(a)*W + (2*sin(a/2)^2)*(W^2);
    rodriguesQuat = rotm2quat(rodriguesRotMatrix);
    rot = quatmultiply(rodriguesQuat, rot);
    
    
    %% Pitch
    u = localPitch; % Set axis for rotation.
    a = pitchAngle; % [rad]
    W = [  0  -u(3)  u(2);
         u(3)   0   -u(1);
        -u(2)  u(1)   0 ];
    I = eye(3); % identity matrix 3x3
    rodriguesRotMatrix = I + sin(a)*W + (2*sin(a/2)^2)*(W^2);
    rodriguesQuat = rotm2quat(rodriguesRotMatrix);
    rot = quatmultiply(rodriguesQuat, rot);
    
    
    %% Yaw
    u = localYaw; % Set axis for rotation.
    a = yawAngle; % [rad]
    W = [  0  -u(3)  u(2);
         u(3)   0   -u(1);
        -u(2)  u(1)   0 ];
    I = eye(3); % identity matrix 3x3
    rodriguesRotMatrix = I + sin(a)*W + (2*sin(a/2)^2)*(W^2);
    rodriguesQuat = rotm2quat(rodriguesRotMatrix);
    rot = quatmultiply(rodriguesQuat, rot);
    
    
    
    
    
    
    
    
    
    
    
%     %% INCLINATION
%     %% Correction for orbit inclination
%     % If the lat-long point is in the ascending portion of the orbit, the 
%     % inclination will be positive.
%     % If the lat-long point is in the descending portion of the orbit, the 
%     % inclination will be negative.
%     
%     % Compute the derivative of the latitude.
%     % If derivative is positive, orbit is ascending.
%     % If derivative is negative, orbit is descending.
%     
%     % Check for positive or negative change in latitude.
%     if i > 1
%       lat_before = parsat(n).lat.signals.values(i-1);
%       lat_now = parsat(n).lat.signals.values(i);
%       diff = lat_now - lat_before;
%     else
%       diff = 0;
%     end
%     
%     % Example A:
%     % If lat_now = 0, lat_before = -1,
%     % diff = +1,
%     % diff positive -> latitude increasing, ascending portion of the orbit
%     % latitude increasing -> use positive orbit inclination for YAW
%     
%     % Example B:
%     % If lat_now = 85, lat_before = 86,
%     % diff = -1,
%     % diff negative -> latitude decreasing, descending portion of the orbit
%     % latitude decreasing -> use negative orbit inclination for YAW
%     
%     if diff < 0
%       inclinationToUse = -parsat(n).inc;
%     else
%       inclinationToUse = parsat(n).inc;
%     end
%     
%     % Set rotation angle on ECEF's X-axis.
%     a = inclinationToUse;
%     rotInc = [cos(a/2), sin(a/2), 0, 0]; % [quaternion]
%     rot = quatmultiply(rotInc, rot);
    
    



%     %% LATITUDE
%     % Set rotation angle on ECEF's Y-axis.
%     a = -lat;
%     rotLat = [cos(a/2), 0, sin(a/2), 0]; % [quaternion]
% %     rot = quatmultiply(rotLat, rot);
%     
%     %% LONGITUDE
%     % Set rotation angle on ECEF's Z-axis.
%     a = long;
%     rotLong = [cos(a/2), 0, 0, sin(a/2)]; % [quaternion]
% %     rot = quatmultiply(rotLong, rot);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     %% ROLL
%     % ROLL on ECEF's Y-axis with angle -latitude.
%     % Output vector2.
% %%    a = -lat + rollAngle;
%     a = -lat;
%     rotRoll = [cos(a/2), 0, sin(a/2), 0]; % [quaternion]
%     
%     % Update locals.
%     localYaw = quatrotate(rotRoll, localYaw);
%     localPitch = quatrotate(rotRoll, localPitch);
%     localRoll = quatrotate(rotRoll, localRoll);
%     
%     % Update cumulative quaternion rotation.
% %%    rot = quatmultiply(rotRoll, rot);
%     
%     %% PITCH
%     % Rotate on ECEF's Z-axis with angle +longitude.
%     % Output vector3.
%     a = long + pitchAngle;
%     rotPitch = [cos(a/2), 0, 0, sin(a/2)]; % [quaternion]
%     
%     % Update locals.
%     localYaw = quatrotate(rotPitch, localYaw);
%     localPitch = quatrotate(rotPitch, localPitch);
%     localRoll = quatrotate(rotPitch, localRoll);
    
    
    
    
    
    % Update cumulative quaternion rotation.
%%    rot = quatmultiply(rotPitch, rot);
    
%     %% Correction for orbit inclination
%     % If the lat-long point is in the ascending portion of the orbit, the 
%     % inclination will be positive.
%     % If the lat-long point is in the descending portion of the orbit, the 
%     % inclination will be negative.
%     
%     % Compute the derivative of the latitude.
%     % If derivative is positive, orbit is ascending.
%     % If derivative is negative, orbit is descending.
%     
%     % Check for positive or negative change in latitude.
%     if i > 1
%       lat_before = parsat(n).lat.signals.values(i-1);
%       lat_now = parsat(n).lat.signals.values(i);
%       diff = lat_now - lat_before;
%     else
%       diff = 0;
%     end
%     
%     % Example A:
%     % If lat_now = 0, lat_before = -1,
%     % diff = +1,
%     % diff positive -> latitude increasing, ascending portion of the orbit
%     % latitude increasing -> use positive orbit inclination for YAW
%     
%     % Example B:
%     % If lat_now = 85, lat_before = 86,
%     % diff = -1,
%     % diff negative -> latitude decreasing, descending portion of the orbit
%     % latitude decreasing -> use negative orbit inclination for YAW
%     
%     if diff < 0
%       inclinationToUse = -parsat(n).inc;
%     else
%       inclinationToUse = parsat(n).inc;
%     end
    
%     %% YAW
%     % Use YAW with Rodrigues rotation on satellite's unit position vector.
%     % ------------------------------------------------------------------------ %
%     % Unit vector from Earth center to satellite location -> u
%     % ------------------------------------------------------------------------ %
%     % References for Rodrigues rotation matrix.
%     % https://math.stackexchange.com/questions/142821
%     % https://mathworld.wolfram.com/RodriguesRotationFormula.html
%     
%     % Set normalized vector 'u'.
%     u = [unitPositionVector2(1) -unitPositionVector2(2) -unitPositionVector2(3)];
%     
%     % Set yaw rotation angle.
%     a = inclinationToUse + yawAngle;
%     
%     % Compute Rodrigues rotation matrix.
%     W = [  0  -u(3)  u(2);
%          u(3)   0   -u(1);
%         -u(2)  u(1)   0 ];
%     I = eye(3); % identity matrix 3x3
%     rodriguesRotMatrix = I + sin(a)*W + (2*sin(a/2)^2)*(W^2);
%     
%     % Convert Rodrigues rotation matrix to quaternions.
%     rotYaw = rotm2quat(rodriguesRotMatrix);
%     
%     % Update locals.
%     localYaw = quatrotate(rotYaw, localYaw);
%     localPitch = quatrotate(rotYaw, localPitch);
%     localRoll = quatrotate(rotYaw, localRoll);
    
    % Update cumulative quaternion rotation.
%%    rot = quatmultiply(rotYaw, rot);
    
    %% Total rotation
    % Multiplied all rotation quaternions to get final rotation.
    parsat(n).rot.signals.values(i,:) = rot;
    
  end
  % Set timer stop time.
  timeProcStop = posixtime(datetime('now')); % Posixtime [seconds].
  timeProcDuration = timeProcStop - timeProcStart;
  msg = sprintf('%3.d/%d ... done. (%.2f seconds)',n,numsats,timeProcDuration);
  send(dq, msg);
  
  % Globally concatenate all output variables on lab index 1.
  % Must be the last lines of code of the parallel pool.
  % satellites = gcat(parsat,1,1);
end % Parallel code.

for i = 1:numsats
  allsats = parsat{i};
  simsat(i) = allsats(i);
end


%% Simulink
% Model names.
fprintf('\n%s','Preparing Simulink models ... ');
% Set timer start time.
timeProcStart = posixtime(datetime('now')); % Posixtime [seconds].
modelMain = 'asbCubeSat';
model3DWorld = 'asbCubeSat/Visualization/Virtual Reality World/';

% Load and open main Simulink model.
open_system(modelMain);


%% Delete all code-generated lines and blocks
blockPaths = find_system(model3DWorld,'Type','Block');

% Find '/_Sat' in blockPaths.
blocksToDelete = find(contains(blockPaths,'/_Sat'));

% Delete lines and blocks.
for n = 1:length(blocksToDelete)
  % Check if block is a main- or sub-block in the Simulink model.
  % First get block name.
  blockName = blockPaths{blocksToDelete(n)};
  % Then check number of '/' occurrences.
  pathSeparators = count(blockName,'/');
  if pathSeparators == 3 % it is a main-block.
    % Delete line(s).
    h = get_param(blockName,'LineHandles');
    for p = 1:length(h.Outport) % for each output port
      if h.Outport(p) > 0 % if line exists
        delete_line(h.Outport(p)); % delete line
      end
    end
    % Delete block.
    delete_block(blockName);
  end
end


%% Update X3D File and VR Block
% Name must match the name of the block in Simulink.
b3DWorld = 'VR Sink1/';

% Open x3D file for writing and discard all contents ('w' option).
file3DWorld = fopen([pathVisualization,filesep,'simulation',filesep,...
  'cosmosSimulation.x3d'],'w');

% Write initial structure before satellites section.
% Open file with initial structure.
file3DWorld_part1 = fopen([pathVisualization,filesep,'simulation',filesep,...
  'cosmosSimulation_struct1.x3d'],'r');
while ~feof(file3DWorld_part1)
  sline = fgetl(file3DWorld_part1);
  fprintf(file3DWorld,'%s\n',sline);
end
% Close file.
fclose(file3DWorld_part1);

% Write N satellite structures.
for n=1:(numsats-1)
  % Define satellite structure.
  structSat = {
    ['  <Transform DEF="Satellite',num2str(n),'"  scale="0.05 0.05 0.05" rotation="0.50838 0.84409 -0.08794 2">']
    ['    <Shape DEF="SatelliteBus',num2str(n),'" >']
    '      <Appearance>'
    '        <ImageTexture url=''"texture/cubesat.jpg" ''>'
    '        </ImageTexture>'
    '        <Material specularColor="1 1 1" shininess="1" emissiveColor="0.2 0.2 0.2" diffuseColor="0.2 0.2 0.2" ambientIntensity="0.6">'
    '        </Material>'
    '      </Appearance>'
    '      <Box size="0.3 0.3 0.3">'
    '      </Box>'
    '    </Shape>'
    ['    <Transform DEF="SolarPanelTwo',num2str(n),'"  translation="0 0.16 -0.3">']
    '      <Shape>'
    '        <Appearance>'
    '          <ImageTexture url=''"texture/solarpanel.jpg" ''>'
    '          </ImageTexture>'
    '          <Material shininess="1" diffuseColor="0.9 0.76733 0.61963" ambientIntensity="0">'
    '          </Material>'
    '        </Appearance>'
    '        <Box size="0.3 0.02 0.3">'
    '        </Box>'
    '      </Shape>'
    '    </Transform>'
    ['    <Transform DEF="SolarPanelOne',num2str(n),'"  translation="0 0.16 0.3">']
    '      <Shape>'
    '        <Appearance>'
    '          <ImageTexture url=''"texture/solarpanel.jpg" ''>'
    '          </ImageTexture>'
    '          <Material shininess="1" diffuseColor="0.9 0.76733 0.61963" ambientIntensity="0">'
    '          </Material>'
    '        </Appearance>'
    '        <Box size="0.3 0.02 0.3">'
    '        </Box>'
    '      </Shape>'
    '    </Transform>'
    ['    <Transform DEF="Antenna',num2str(n),'" >']
    ['      <Transform DEF="AntennaDish',num2str(n),'"  translation="0 -0.15 0">']
    '        <Shape>'
    '          <Appearance>'
    '            <Material specularColor="1 1 1" shininess="1" diffuseColor="0.9 0.65504 0.31917">'
    '            </Material>'
    '          </Appearance>'
    '          <Cone height="0.07" bottomRadius="0.06">'
    '          </Cone>'
    '        </Shape>'
    '      </Transform>'
    '    </Transform>'
    ['    <Transform DEF="FollowSat',num2str(n),'"  translation="-5.5 2 0" rotation="-0.0923 -0.7011 -0.0923 1.5878">']
    ['      <Viewpoint DEF="FollowSatView',num2str(n),'"  fieldOfView="0.7854" description="Lock Satellite-',num2str(n),'" position="0 0 0">']
    '      </Viewpoint>'
    '    </Transform>'
    '  </Transform>'
  };
  
  % Write structure to X3D file.
  fprintf(file3DWorld,'%s\n',structSat{:});
end

% Define final structure after satellites section.
structFinal = {
  '</Transform>'
  '</Scene>'
  '</X3D>'
};

% Write structure to X3D file.
fprintf(file3DWorld,'%s\n',structFinal{:});

% Close x3D file.
fclose(file3DWorld);

% Help!
% VR block input ports - To see structure of the input ports:
% get_param([model3DWorld,b3DWorld],'FieldsWritten')





%%TO DO
% Define length of variable holding the description of the input ports.





% Define input ports order.
b3DWorldInputPorts = [...
  'Earth.rotation.4.1.1.double',...
  '#Sun.translation.3.1.1.double',...
  '#Umbra.rotation.4.1.1.double',...
  '#Umbra.translation.3.1.1.double',...
  '#FollowReference.rotation.4.1.1.double',...
  '#FollowReference.translation.3.1.1.double',...
  '#SatelliteReference.rotation.4.1.1.double',...
  '#SatelliteReference.translation.3.1.1.double'];
for n=1:(numsats-1)
  b3DWorldInputPorts = [b3DWorldInputPorts,...
    '#Satellite',num2str(n),'.rotation.4.1.1.double',...
    '#Satellite',num2str(n),'.translation.3.1.1.double'];
end

% Update input ports to VR block.
set_param([model3DWorld,b3DWorld],'FieldsWritten',b3DWorldInputPorts);


%% Add Blocks and Lines
% Update input variables for reference satellite.
refSatVar = ['simsat(',num2str(numsats),')'];
set_param([model3DWorld,'Local_Ref_Rot'],'VariableName',[refSatVar,'.rot']);
set_param([model3DWorld,'Local_Ref_Pos'],'VariableName',[refSatVar,'.pos']);

% Add blocks and lines for each satellite in the formation.
for n=1:(numsats-1)
  % Help!
  % To check the dialog parameters of a block:
  % get_param(blockPath,'DialogParameters')
  bSatPos = ['_Sat',num2str(n),'_Pos'];
  add_block('simulink/Sources/From Workspace',[model3DWorld,bSatPos],...
    'VariableName',['simsat(',num2str(n),').pos']);
  
  bSatPosConv = ['_Sat',num2str(n),'_ECEFtoVRML'];
  add_block([model3DWorld,'ECEF to VRML'],[model3DWorld,bSatPosConv]);
  
  bSatRot = ['_Sat',num2str(n),'_Rot'];
  add_block('simulink/Sources/From Workspace',[model3DWorld,bSatRot],...
    'VariableName',['simsat(',num2str(n),').rot']);
  
  bSatRotConv = ['_Sat',num2str(n),'_QUATtoVRML'];
  add_block([model3DWorld,'Quaternion to VRML'],[model3DWorld,bSatRotConv]);
  
  % Help!
  % To get the port connectivity of the 3D World block:
  % p = get_param([model3DWorld,b3DWorld],'PortConnectivity')
  % for n=1:length(p), disp([num2str(n),' : ',num2str(p(n).Position)]), end
  
  % Compute 3D World's input port numbers for current satellite.
  pSatN = 9 + (n-1)*2;
  pRotation    = num2str(pSatN);
  pTranslation = num2str(pSatN + 1);
  
  
  
  
  
  
  
  %PUT scale directly into xml code generation for satellites
  satSizeInitial = [2e5 2e5 2e5]; % defined by Simulink model.
  satSizeConverted = satSizeInitial .* (1/.3e6); % Model to VRML conversion.
  satUserScale = 0.1; % defined by user (configurable).
  satSizeFinal = satSizeConverted .* satUserScale; % [1x3 VRML]
  % 0.066666666666667 -> 0.067
    
  
  
  
  
  
  
  
  % Add lines to connect blocks.
  add_line(model3DWorld,[bSatPos,'/1'],[bSatPosConv,'/1']);
  add_line(model3DWorld,[bSatPosConv,'/1'],[b3DWorld,pTranslation]);
  
  add_line(model3DWorld,[bSatRot,'/1'],[bSatRotConv,'/1']);
  add_line(model3DWorld,[bSatRotConv,'/1'],[b3DWorld,pRotation]);
end

% Auto-arrange model3DWorld.
% Not necessary yet.

% Save model.
save_system(modelMain);
% Set timer stop time.
timeProcStop = posixtime(datetime('now')); % Posixtime [seconds].
timeProcDuration = timeProcStop - timeProcStart;
fprintf('done. (%.2f seconds)\n',timeProcDuration);


%% Check and Run
% Help!
% To see all object handles open in MATLAB, enter in Command Window:
% object_handles = findall(groot)
% To find handle of the figure for 3D visualization, enter in Command Window:
% findall(groot,'Name','COSMOS Visualization')

% Create a virtual world associated with the virtual world 3D file filename and
% return its handle. If the VR world already exists, a handle to the existing
% virtual world is returned. Specify the file name as a string.
vrWorldObj = vrworld('cosmosSimulation.x3d');
reload(vrWorldObj);

% Get handle of the current virtual reality figure. The current virtual reality
% figure is the currently active virtual reality figure window in which you can
% get and set the viewer properties. If no virtual reality figure exists, the
% MATLAB software returns an empty vrfigure object.
% vrFigHandle = vrgcf;

% Try to get handle of the figure window for 3D visualization.
vrFigHandle = findall(groot,'Name','COSMOS Visualization');
if isempty(vrFigHandle)
  % If handle is empty, means that the virtual world figure is closed.
  % Then, open VR figure.
  vrFigHandle = vrfigure(vrWorldObj);
  % cosmosVisHandle = findall(groot,'Name','COSMOS Visualization');
  % Find VRSimMenu in cosmosVisHandle.Children(5)
  % Find submenus in cosmosVisHandle.Children(5).Children(1) to (3)
else
  % If handle is found, bring figure to front.
  figure(vrFigHandle);
end

% Set timer stop time.
timeScriptStop = posixtime(datetime('now')); % Posixtime [seconds].
timeScriptDuration = timeScriptStop - timeScriptStart;
fprintf('\nTotal visualization processing time: %.2f seconds\n',timeScriptDuration);

fprintf('%s','Ready to play visualization: ')
if AUTORUN
  fprintf('%s\n','Autoplay [ON]')
  % Start Simulink visualization.
  sim(modelMain, stopTime);
else
  fprintf('%s\n','Autoplay [OFF]')
  % Confirmation here:
  msgfig = msgbox('Visualization Processing Completed','MATLAB Info','help','modal');
  uiwait(msgfig);
end

% Clear composites.
% Needed to close active spmd.
% See: https://nl.mathworks.com/matlabcentral/answers/223939-how-to-specify-the-number-of-labs-in-smpd
clear a adeg aEarthDeg aEarthRad angPos angX c I lat localPitch localRoll localYaw long msg ...
  normalOrbitVector normalVector pitchAngle pitchCorrection rodriguesQuat ...
  rodriguesRotMatrix rollAngle rot rotEarthFix rotInc ...
  sma time timeProcDuration timeProcStart timeProcStop u v1 v2 W x xyzUnit ...
  y yawAngle z satellites parsat;
