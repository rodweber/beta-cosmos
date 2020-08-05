%% File to run Cosmos software

%% Set paths and MATLAB parameters

warning on verbose;
delete(gcp('nocreate'));
close all; clear all; clc; %#ok<CLALL>

% Inform the name of this file without the extension "m".
THIS_FILE_NAME = 'runCosmosBeta';

if(~isdeployed)
	
	% Get directory path of the active file in MATLAB's Editor.
	[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
	addpath(filepath); % Add file path to the current MATLAB path.
	
	% Get directory path of the file name set in THIS_FILE_NAME.
	[filepath,~,~] = fileparts(which(THIS_FILE_NAME));
	addpath(filepath); % Add file path to the current MATLAB path.
	
	% Change working directory to the directory of this m-file.
	cd(filepath);
	
end

%% Set parameters for the simulation

% - Struct of parameters:
%   NumSatellites    : Total number of satellites in the formation.
%   FormationMode    : Mode for the satellites formation flight.
%   Altitude         : Height above sea level [m].
%   DeltaAngle       : Roll, pitch, yaw angles resolution [deg].
%   AutoResponse     : If satellite should send responses [bool].
%   AvailableGPS     : GPS availability [bool].
%   AvailableTLE     : TLE availability [bool].
%   MaxNumOrbits     : Maximum number of orbits to run.
%   OrbitSectionSize : Size of each orbit section [deg].
%   InitIDX          : Initial idx.
%   AccelFactor      : Acceleration factor for the simulation.

parametersIvanov = struct( ...
	'NumSatellites'   , 4     , ...
	'FormationMode'   , 1     , ...
	'Altitude'        , 340e3 , ...
	'DeltaAngle'      , 30    , ...
	'AutoResponse'    , true  , ...
	'AvailableGPS'    , true  , ...
	'AvailableTLE'    , false , ...
	'MaxNumOrbits'    , 2     , ...
	'OrbitSectionSize', 2     , ...
	'InitIDX'         , 1     , ...
	'AccelFactor'     , 100e3 );

% where are?:
%satellite mass
%number of solar panels and their size etc
%eclipsedSun=0;
%rotatingSun=1;
parametersCLUSTER = struct( ...
	'NumSatellites'   , 3     , ...
	'FormationMode'   , 1     , ...
	'Altitude'        , 453e3 , ...
	'DeltaAngle'      , 45    , ...
	'AutoResponse'    , true  , ...
	'AvailableGPS'    , true  , ...
	'AvailableTLE'    , false , ...
	'MaxNumOrbits'    , 5     , ...
	'OrbitSectionSize', 0.2   , ...
	'InitIDX'         , 1     , ...
	'AccelFactor'     , 100e3 , ...
  'vizScale'        , 1     , ...
  'FolderFFPS'      , 'ffps', ...
  'FFPS'            , [struct('ffp1',0,...
                              'ffp2',0,...
                              'ffp3',0,...
                              'ffp4',0,...
                              'ffp5',0,...
                              'ffp6',0,...
                              'ffp7',0,...
                              'ffp8',0)
                       struct('ffp1',0,...%% normally 100
                              'ffp2',-115,...
                              'ffp3',0,...%% normally 100
                              'ffp4',0,...
                              'ffp5',0,...
                              'ffp6',0,...
                              'ffp7',0,...
                              'ffp8',0)
                       struct('ffp1',0,...
                              'ffp2',-2*115,...
                              'ffp3',0,...
                              'ffp4',0,...
                              'ffp5',0,...
                              'ffp6',0,...                              
                              'ffp7',0,...
                              'ffp8',0)      
                             ]);
parametersISMission = struct( ...
	'NumSatellites'   , 2     , ...
	'FormationMode'   , 1     , ...
	'Altitude'        , 453e3 , ...
	'DeltaAngle'      , 45    , ...
	'AutoResponse'    , true  , ...
	'AvailableGPS'    , true  , ...
	'AvailableTLE'    , false , ...
	'MaxNumOrbits'    , 1     , ...
	'OrbitSectionSize', 0.1   , ...
	'InitIDX'         , 0     , ...
	'AccelFactor'     , 100e3 );

iniConditionsIvanov=[0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0];
iniConditionsCLUSTER=[0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0];
iniConditionsISMission=[0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0];

%% this.SolarFactor,this.WindFactor, this.SatelliteMass ,this.Panels,this.SurfacePanel,this.SurfaceRef can be found in FlightControl
%% desired values are in FlightControl
%% initial conditions are in Satellite.initialize
%% parameters for riccati equation can be found in riccatiequation

%% Instantiate a simulation object with the selected parameters.
% Don't use variable name 'sim'. It is a reserved name needed for Simulink.
%csim = CosmosSimulation(parametersIvanov,iniConditionsIvanov);
csim = CosmosSimulation(parametersCLUSTER,iniConditionsCLUSTER);
%csim = CosmosSimulation(parametersISMission,iniConditionsISMission);

%% Initiate and run simulation.
csim.startSimulation();

%% create objects needed for documentation generation
% Create global alias for the array of satellites.
sat = csim.Satellites; % Aliases: sat(1) to sat(n).
% Create global alias for the array of orbits.
orbit = csim.Orbits; % Aliases: orbit(1) to orbit(n).
% Create global alias for the array of flight control modules.
fc = csim.FlightControlModules; % Aliases: fc(1) to fc(n).
% Create global alias for the array of GPS modules.
gps = csim.GPSModules; % Aliases: gps(1) to gps(n).



%% Plot and visualize results
%! this should go to sim.start
%! empty place holders for now
csim.plotting(csim.NumSatellites, orbit(1).MeanMotionRad);
csim.visualizationLONLATALT(csim.NumSatellites,orbit(1).Altitude)
%Before:
%u = zeros(csim.NumSatellites, 3, size(csim.SatStates,3));
%e = zeros(csim.NumSatellites, 6, size(csim.SatStates,3));
%csim.plotting(csim.SatStates(:,7:9,:), csim.SatStates(:,1:6,:), csim.SatPositions(1,:,:), csim.TimeVector', csim.NumSatellites, orbit(1).MeanMotionRad, u, e);
%csim.visualizationLONLATALT(csim.NumSatellites,csim.TimeVector',squeeze(csim.SatStates(:,1,:)),squeeze(csim.SatStates(:,2,:)),csim.SatStates(:,3,:),squeeze(csim.SatStates(:,7,:)),squeeze(csim.SatStates(:,8,:)),squeeze(csim.SatStates(:,9,:)),orbit(1).Altitude)
%n visualizationLONLATALT(ns                ,  cosmostime   ,sstx                         ,ssty                         ,sstz                ,pitch                        ,yaw                          ,roll                         ,altitude)



%% save data for documentation generation
% Save current MATLAB workspace variables.
warning off parallel:lang:spmd:CompositeSave;
workspaceFileName = 'workspace.mat';
save(fullfile(filepath, workspaceFileName));
% Print custom objects and classes used.
csim.createListCustomClasses(filepath, workspaceFileName);

% For reference, this was before:
%%% Plot results
% angles = csim.SatStates(:,7:9,:);
% sst = csim.SatStates(:,1:6,:);
% refPosChange = csim.SatPositions(1,:,:);
% time = csim.TimeVector';
% ns = csim.NumSatellites;
% meanMotion = orbit.MeanMotionRad;
% u = zeros(csim.NumSatellites, 3, size(csim.SatStates,3));
% e = zeros(csim.NumSatellites, 6, size(csim.SatStates,3));
% csim.plotting(angles, sst, refPosChange, time, ns, meanMotion, u, e);

fprintf('\nDone.\n\n');