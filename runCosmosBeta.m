%% File to run Cosmos software

%% Set paths and MATLAB parameters
clc;
warning on verbose;
delete(gcp('nocreate'));

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
	
	% Change working directory to the directory of this m-file.
	cd(filepath);
end

%% Instantiate a simulation object
% Do not use 'sim' as variable name. It is a reserved name in Simulink.
csim = CosmosSimulation();

%% Start simulation
csim.startSimulation();

%% Confirmation
fprintf('\nDone.\n\n');
msgfig = msgbox('CosmosBeta Simulation Completed','MATLAB Info','help','modal');
uiwait(msgfig);
