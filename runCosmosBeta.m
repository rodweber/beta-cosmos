%% File to run Cosmos software
% Use code line below only for performance tests.
% proj = currentProject; close(proj); close all; clear all;

% JT Notes
% this.SolarFactor,this.WindFactor, this.SatelliteMass ,this.Panels,
%   this.SurfacePanel,this.SurfaceRef can be found in FlightControl
% desired values are in FlightControl
% initial conditions are in Satellite.initialize
% parameters for riccati equation can be found in riccatiequation

%% Set paths and MATLAB parameters
clc;
warning on verbose;
% Try to get handle of the figure window for 3D visualization.
cosmosVisHandle = findall(groot,'Name','COSMOS Visualization');
% If handle for 3D visualization is empty, simply close all figures.
if isempty(cosmosVisHandle)
  close all;
else % If handle is found, close all figures with exception of 3D visualization.
  cab(cosmosVisHandle);
end

% Set timer start time.
disp('Preparing MATLAB environment ...');
timeMatlabStart = posixtime(datetime('now')); % Posixtime [seconds].

% Inform the name of this file without the extension "m".
THIS_FILE_NAME = 'runCosmosBeta';

if(~isdeployed)
  % Get directory path of the active file in MATLAB's Editor.
  [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
  addpath(filepath); % Add file path to the current MATLAB path.

  % Get directory path of the file name set in THIS_FILE_NAME.
  [filepath,~,~] = fileparts(which(THIS_FILE_NAME));
  addpath(filepath); % Add file path to the current MATLAB path.

  % Add path to ancillary folders.
  addpath(strcat(filepath,filesep,'config'));
  addpath(strcat(filepath,filesep,'utils'));
  addpath(strcat(filepath,filesep,'Experiment-Times'));

  % Change working directory to the directory of this m-file.
  cd(filepath);
end

%% Read parameters from configuration files
% Read and set initial condition parameters from configuration files.
% To change a parameter's value, open the respective JSON file and modify it.
configFF  = readjson('configFormationFlight.json');
configSim = readjson('configSimulation.json');
configVis = readjson('configVisualization.json');

% Get identifier of the selected formation flight mode.
selectedModeID = configFF.SelectedFormationFlightMode;

% Go through list of formation flight modes.
for i = 1:length(configFF.Modes)
  if configFF.Modes(i).ModeID == selectedModeID
    selectedMode = configFF.Modes(i);
    break;
  else
    %% Throw error that Mode ID could not be found.
  end
end

%% Instantiate a simulation object
% Do not use 'sim' as variable name. It is a reserved name in Simulink.
csim = CosmosSimulation(...
  configFF.NumberOfOrbitsForMode,...
  configSim.SimAccelerationFactor,...
  selectedMode.ModeID,...
  selectedMode.ModeName,...
  selectedMode.NumSatellites,...
  selectedMode.SatelliteMass,...
  selectedMode.PanelArea,...
  selectedMode.PanelQuantity,...
  selectedMode.AttitudeResolutionDeg,...
  selectedMode.OrbitInitAltitude,...
  selectedMode.OrbitSectionSizeDeg,...
  selectedMode.OrbitInitSectionID,...
  selectedMode.RMatrixDiagonal,...
  selectedMode.AutoResponse,...
  selectedMode.AvailableGPS,...
  selectedMode.AvailableTLE,...
  selectedMode.InitialConditions,...
  selectedMode.FFPSFolder,...
  selectedMode.FFPSValues);

%% Start simulation
csim.startSimulation();

%% Set global object aliases for autonomous documentation generation
sat   = csim.Satellites; % Aliases: sat(1) to sat(n).
orbit = csim.Orbits; % Aliases: orbit(1) to orbit(n).
fc    = csim.FlightControlModules; % Aliases: fc(1) to fc(n).
nav   = csim.GPSModules; % Aliases: nav(1) to nav(n).
comms = csim.CommModules; % Aliases: comms(1) to comms(n).

%% Post-processing and plotting methods
% Use a Kepler propagator to compute ECEF coordinates

%!RW: Move this to openvis
csim.ECEFprocessing(...
  selectedMode.OrbitInitAltitude,...
  selectedMode.OrbitInclinationDeg,...
  selectedMode.OrbitLANDeg,...
  configVis.VizScale,...
  configVis.VizTimeStep,...
  configVis.V0,...
  false,...
  true,...
  configSim.ParentFolderForCoordFolders,...
  configSim.FolderForNonScaledLLRRPY,...
  configSim.FolderForScaledLLRRPY,...
  configSim.FolderForScaledXYZRPY,...
  configSim.EnablePlotting);
%^^

% Compute the location of the specular point for the GNSS-reflectometry method
if(configSim.EnableGNSSRProcessing)
  csim.GNSSRprocessing(...
    csim.NumSatellites,...
    orbit(1).MeanEarthRadius,...
    configSim.ParentFolderForCoordFolders,...
    configSim.FolderForNonScaledLLRRPY);
end

% Plot and visualize results
if(configSim.EnablePlotting)
  csim.plotting(csim.NumSatellites, orbit(1).MeanMotionRad);
end

%% Save data for autonomous documentation generation
% Save current MATLAB workspace variables.
warning off parallel:lang:spmd:CompositeSave;
% Set name for file with workspace variables.
workspaceFileName = '.workspace.mat';
% Get a list of all variables.
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object.
tosave = cellfun( @isempty, regexp(...
  {allvars.class},...
  '^(matlab\.(ui|graphics|project)|Simulink|parallel|vrworld)') );
% Pass these variable names to save.
save(fullfile(filepath, workspaceFileName), allvars(tosave).name);
% Print custom objects and classes used.
csim.createListCustomClasses(filepath, workspaceFileName);

%% Autorun 3D visualization
if(configSim.AutoOpen3DVisualization)
  % Confirmation inside openvis:
  run('openvis.m');
end

%% End
% Set timer stop time.
timeMatlabStop = posixtime(datetime('now')); % Posixtime [seconds].
timeMatlabDuration = timeMatlabStop - timeMatlabStart;
fprintf('\nDone. Runtime: %.3f seconds\n\n',timeMatlabDuration);

% Automatically open 3D visualization.
if ~configSim.AutoOpen3DVisualization
  % Confirmation here:
  msgfig = msgbox('CosmosBeta Simulation Completed','MATLAB Info','help','modal');
  uiwait(msgfig);
end
