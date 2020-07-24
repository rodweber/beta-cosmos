%% File to open visualization for Cosmos Beta in MATLAB Simulink

%% Set parameters
warning on verbose;
delete(gcp('nocreate'));
close all; clc; % Do not clear variables, it makes the project unstable.

% Set parameter to automatically run Simulink visualization.
AUTORUN = 1; % [true: 1 | false: 0]

% Inform the name of this file without the extension "m".
THIS_FILE_NAME = 'openvis';

% Inform the name of the project file for the Cosmos simulation.
PROJECT_FILE_NAME = 'cosmosVisualization'; % Without extension 'prj'.
PROJECT_FOLDER = 'visualization';

% Files with coordinates data for satellites.
coordfiles = {'sat1LLR'};

%% Set paths
if(~isdeployed)
	% Get directory path of the active file in MATLAB's Editor.
	[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
	addpath(genpath(filepath)); % Add file path to the current MATLAB path.
	
	% Get path for directory of the file name set in THIS_FILE_NAME.
	[filepath,~,~] = fileparts(which(THIS_FILE_NAME));
	addpath(genpath(filepath)); % Add file path to the current MATLAB path.
    
    % Get path for directory of the simulation project.
    [projectPath,~,~] = fileparts(which([PROJECT_FILE_NAME,'.prj']));
	
	% Change working directory to the directory of this m-file.
	cd(filepath);
end

%% Run Simulink file

% Check if there is any Matlab project already open.
if(isempty(matlab.project.rootProject))
    % If there is no project open, launch project from path.
    proj = openProject(projectPath);
else
    % If there is a project already open, get the Project object.
    proj = currentProject();
    % Check if the name of the project is correct.
    if(~strcmp(proj.Name,PROJECT_FILE_NAME))
        % If name is different, show error.
        error(['Project ''%s'' is already open. ',...
            'Close it before running ''%s''.'], ...
            proj.Name, PROJECT_FILE_NAME);
    else
        % If name is correct, show message.
        fprintf(2,'Project with name ''%s'' already open.\n',proj.Name);
        fprintf('To prevent issues, check if project path is correct.\n');
        fprintf('Project path: %s\n',proj.RootFolder);
    end
end

% Define satellite struct.
simsat(1).lat = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).long = struct('time',0,'signals',struct('dimensions',0,'values',0));
simsat(1).sma = struct('time',0,'signals',struct('dimensions',0,'values',0));

% Allocate memory for array of satellites.
numsats = length(coordfiles);
simsat = repmat(simsat(1),numsats,1);

for n = 1:numsats
    coord = readmatrix(coordfiles{n});
    timestamps = coord(:,1); % [seconds]
    latDeg = coord(:,2); % [degrees]
    longDeg = coord(:,3); % [degrees]
    smaKm = coord(:,4); % [km]
    
    % Convert data to SI units
    latRad = latDeg * (pi/180); % [rad]
    longRad = longDeg * (pi/180); % [rad]
    smaMeters = smaKm * 1000; % [m]
    
    % Place data into satellite struct.
    simsat(n).lat.time = timestamps;
    simsat(n).lat.signals.dimensions = 1;
    simsat(n).lat.signals.values = latRad;
    
    simsat(n).long.time = timestamps;
    simsat(n).long.signals.dimensions = 1;
    simsat(n).long.signals.values = longRad;
    
    simsat(n).sma.time = timestamps;
    simsat(n).sma.signals.dimensions = 1;
    simsat(n).sma.signals.values = smaMeters;
end

% Get last timestamp of the time vector, set it as simulation time.
stopTime = simsat(1).lat.time(end);

% Check length of lat and long vectors.

% Get each pair of lat-long.

    % Set original vector1 as [SMA, 0, 0].
    
    % Correct orientation of the 3D satellite.
    % Rotate -90deg (-pi/2) on ECEF's Y-axis.
    
    % Rotate on ECEF's Y-axis with angle -latitude.
    % Output vector2.
    
    % Rotate on ECEF's Z-axis with angle +longitude.
    % Output vector3.
    
    
    
    
    
    
    
    % ---------------------------------------------------------------------
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
    % ---------------------------------------------------------------------
    

    % Method 2
    % For each pair of lat-long, compute vector to center of Earth.
    % Start with vector1 = [1 0 0].
    
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
    
    
    
% Change working directory to the directory of the project file.
cd(fullfile(filepath,PROJECT_FOLDER));

% Open Simulink model and run it.
open_system('asbCubeSat');

% To see all object handles open in MATLAB, enter in Command Window:
% object_handles = findall(groot)
% To find object handle of the figure window for 3D visualization, enter
% in Command Window:
% findall(groot, 'Name', 'COSMOS Visualization')

if AUTORUN
    % Try to get handle of the figure window for 3D visualization.
    cosmosVisHandle = findall(groot, 'Name', 'COSMOS Visualization');
    if isempty(cosmosVisHandle)
        % If handle is empty, means that the figure is closed.
        % Open figure.
        open('cosmosSimulation.x3d');
        cosmosVisHandle = findall(groot, 'Name', 'COSMOS Visualization');
        % Find VRSimMenu in cosmosVisHandle.Children(5)
        % Find submenus in cosmosVisHandle.Children(5).Children(1) to (3)
    else
        % If handle is found, figure is already open. Bring it to front.
        figure(cosmosVisHandle);
    end
    % Start Simulink visualization.
    sim('asbCubeSat',stopTime);
end
