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

%% JT Notes
% this.SolarFactor,this.WindFactor, this.SatelliteMass ,this.Panels,
%   this.SurfacePanel,this.SurfaceRef can be found in FlightControl
% desired values are in FlightControl
% initial conditions are in Satellite.initialize
% parameters for riccati equation can be found in riccatiequation


%% Instantiate a simulation object
% Do not use 'sim' as variable name. It is a reserved name in Simulink.
csim = CosmosSimulation();

%% Start simulation
csim.startSimulation();

%% Set global object aliases for autonomous documentation generation
sat   = csim.Satellites; % Aliases: sat(1) to sat(n).
orbit = csim.Orbits; % Aliases: orbit(1) to orbit(n).
fc    = csim.FlightControlModules; % Aliases: fc(1) to fc(n).
gps   = csim.GPSModules; % Aliases: gps(1) to gps(n).


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


%% Save data for autonomous documentation generation
% Save current MATLAB workspace variables.
warning off parallel:lang:spmd:CompositeSave;
workspaceFileName = 'workspace.mat';
save(fullfile(filepath, workspaceFileName));
% Print custom objects and classes used.
csim.createListCustomClasses(filepath, workspaceFileName);

%!RW: For reference, this was before:
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
