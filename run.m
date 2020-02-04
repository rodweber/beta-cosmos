%% File to run cosmosFS.
%
% Priority:
% - @cosmos & @cosmosFS: first uses riccati outside the loop, another
%   uses riccati inside the loop
% - @cosmosFS.m,
%   line 241: orbitalproperties outputs mean motion in [rad/s]; 
%   line 243: whereInWhatOrbit converts to [deg/s];
%   line 95: converts back to [rad/s]
%   @riccatiequation, line 7: mean motion should be in deg/s ?
% - @cosmosFS.m, lines 39, 92: what is idx? idx in @whereInWhatOrbit?
%   is it the ID for the section? why starts at 120?
% - @cosmosFS.m, line 234: check possible values for if, does it
%   make sense? does it always go if = true?
% - go over all steps from MYcosmosFS.m (spmd loop): now at line 136
%
% To do:
% - Add docs('update') option to update publish for all m files
% - Generate pdf publish files for all files in Windows PC
% - Upload to Git all publish files generated in Windows PC
% - Add docs('all') to open publish files as well as docs for classes
% - Update readme.md
%   - To run CosmosFS program, open and run the file 'run.m' in MATLAB
%   - To see documentation for the custom objects and classes used, 
%     enter command 'docs' in MATLAB's Command Window
%   - To see documentation for all files in CosmosFS program, 
%     enter command 'docs('all')' in MATLAB's Command Window
%   - To generate PlantUML code for CosmosFS class diagram, 
%     enter command 'uml' in MATLAB's Command Window
% - check usage of var wind
% - check usage of var refSurf
% - review @aeroPressureForce.m
% - review @aeroDragLiftSentman.m
% - review @vectorRotation.m
% - review @solarPressureForce.m
% - Use class property attribute 'Constant, GetAccess = public' ...
%   for constant attributes e.g. in class Orbit
% - Change class IvanovFormationFlight to FormationFlight, add
%   property 'case' and set property as 'Ivanov'
% - Change class IvanovSatellite to Satellite, add
%   property 'case' and set property as 'Ivanov'
% - @orbitalproperties.m, line 89: Check function semi-major axis 
%   for possible simplification
% - Create package '+cosmos' for all custom functions
% - Custom package to be downloaded and saved into users MATLAB folder
% - Create tutorial to place custom package into users MATLAB folder
%   for both Windows and Mac
%
% Recently done:
% - [3] Add orbit section loop into parloop
% - [2] Add idx and pause into parloop
% - [1] Change calls for function path to addpath
% - Fix order of commands in function uml (2)
% - ^3 Fix list of git commands to update branch 'out'
% - ^2 Update function uml for Windows compatibility (2)
% - Add branch 'dev' and set uml output to branch 'out'
% - Use proxy service of the PlantUML server to open UML diagram
% - Add temp.uml to GitHub automatically (3)
% - Remove hyperlinks, tooltips and fix visibility codes in m2uml
% - Add spacing between VisibilityCode and Names in m2uml
% - In +m2uml.Property, lines 100 to 109: remove Hyperlink, ToolTip
% - In +m2uml.Operation, lines 87 to 94: remove Hyperlink, ToolTip
% - In +m2uml.ClassNode, lines 134 to 144: remove Hyperlink, ToolTip
% - Fix filesep in m2uml.filespec2fqn for dos/unix compatibility (2)
% - Update custom function uml
% - Fix function publish in function docs as an optional parameter
% - Fix function docs to work in both Windows and Mac
% - Fix documentation tool that shows custom object classes used
% - Fix change of working directory and path of the running m-file
% - Change method to get full path of file run
% - Add class Main and change main.m to run.m
% - Check MeanAnomalyFromAN value
% - Add riccati equation
% - Add to class Orbit the property mean anomaly from ascending node
% - Add to class Orbit the property of mean motion in rad/s
% - Fix implementation of whereInWhatOrbit
% - Add getCurrentOrbitNumber to class Satellite
% - Remove redundant properties in class IvanovFormationFlight
% - Check orbital properties
% - Add function to update orbital parameters in class Orbit
% - Add variables to set GPS/TLE availability
% - Fix orbit counter increment and checkpoint
% - Add print checkpoints in main and end of files
% - Add getStatus() to class Simulation
% - Add orbit counter to IvanovSatellite
% - Add set and get methods for formation flight mode
% - Add class Orbit into IvanovFormationFlight
% - Add class CosmosSimulation to main file
% - Add lib folder under the same directory of the main file
% - Add code to automatically update the working directory
% - Add array of IvanovSatellite objects into IvanovFormationFlight


%% Parameters

warning on verbose;
close all; clear all; clc; %#ok<CLALL>

% Inform the name of this file without the extension "m".
THIS_FILE_NAME = 'run';

if(~isdeployed)
	
	% Get directory path of the active file in MATLAB's Editor.
	[filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
	addpath(filepath); % Add file path to the current MATLAB path.
	
	% Get directory path of the file name set in THIS_FILE_NAME.
	[filepath,~,~] = fileparts(which(THIS_FILE_NAME));
	addpath(filepath); % Add file path to the current MATLAB path.
	
	% Change working directory to the directory of this m-file.
	cd(filepath);
	
	% Add lib folder to the current MATLAB path.
	current_path = path;
	path(current_path,strcat('.',filesep,'lib',filesep));
	
end

% Instantiate object of class CosmosSimulation.
max_number_of_orbits = 2;
acceleration_factor = 10000;
available_GPS = false;
available_TLE = false;
sim = CosmosSimulation(max_number_of_orbits,acceleration_factor);

% Instantiate object of class Orbit.
altitude = 340000; % [meters].
orbit = Orbit(altitude);

% Instantiate object of class IvanovFormationFlight.
number_of_satellites = 4;
iv = IvanovFormationFlight(orbit,number_of_satellites,...
                           available_GPS,available_TLE);

% Create aliases for satellite objects.
sat = iv.Satellites; % Aliases: sat(1) to sat(n).

% Initial idx.
idx = 120;
orbitSectionSize = 2; % Size of each orbit section [deg].
orbitSections = 1:orbitSectionSize:360;


%% Parallel loop

% Create data queue for parallel pool.
dq = parallel.pool.DataQueue;

% Define function to call when new data is received on the DataQueue.
afterEach(dq, @disp);

% Create and return pool with the specified number of workers.
parpool(number_of_satellites);

% Set the start time for the parallel pool.
startTimePool = posixtime(datetime('now')); % Posixtime [seconds].

% Execute parallel code on workers of parallel pool.
spmd(number_of_satellites)
	
	% Set satellite IDs (id) for each of the satellites.
	id = labindex;
	alive = true;
	send(dq,['Satellite number ',num2str(id),' is alive.']);
	
	% Sattelites are alive but still doing nothing.
	while alive
		
		% Get formation flight mode.
		fofliMode = iv.getFormationFlightMode();
		
		% Get current orbit number of the satellite.
		currentOrbit = sat(id).getCurrentOrbitNumber();
		
		% Get updated simulation status: 0 = Stop; 1 = Good.
		[sim_status, msg] = sim.getStatus(currentOrbit);
		
		% Log.
		% send(dq,['Sat.',num2str(id),': orbit counter = ',...
		% num2str(currentOrbit),'. ',msg]);
		
		% Get battery status from the satellite.
		battery_status = 1;
		if battery_status
			% Switch on the GPS.
		end
		
		% Orbital loop.
		while sim_status % While simulation status is all good.
			
			% Increment the orbit counter of the satellites.
			sat(id).incrementOrbitCounter();
			
			% Set the start time for the current satellite orbit.
			startTimeOrbit = posixtime(datetime('now')); % Posixtime [s].
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% Calculate endOfSectionsCycle.
			endOfSectionsCycle = (idx-1)/size(orbitSections,2);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% Update orbital parameters for the satellites.
			sat(id).whereInWhatOrbit(endOfSectionsCycle);
			
			% Log.
			send(dq,['Sat ',num2str(id),' - MeanAnomalyFromAN = ',...
			         num2str(sat(id).Orbit.MeanAnomalyFromAN)]);
			
			% Settings for control algorithm, is this necessary every orbit?
			[P,IR,A,B] = riccatiequation(orbit.MeanMotionRad,iv.SSCoeff);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% Wait until end of orbit section.
			idx = find(orbitSections >= orbit.MeanAnomalyFromAN,1,'first');
			idx = idx + 1;
			
			pause((orbitSections(idx) - orbit.MeanAnomalyFromAN) / ...
				orbit.MeanMotion / sim.AccelFactor);
			
			% Orbit sections loop.
			while sim_status && (idx <= size(orbitSections,2))
				
				% Determine cycle start time in order to allow subtraction of 
				% the cycle duration from waiting.
				%startTimeSection = posixtime(datetime('now')); % Posixtime [s].
				%startTimeSection = now();
				
				% To do:
				% Set attitude computed in last iteration.
				
				% To do:
				% Compute attitude for next section.
				
				% Determine desired trajectory.
				
				
				
				
				
				
				
				
				
				
				
				
				% Increment section counter.
				idx = idx + 1;
				
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%% MORE CODE HERE %%%%%%%%%%%%%%%%%%%%%%%
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				%
				
			end % Orbital sections loop.
			
			% Get current orbit number of the satellite.
			currentOrbit = sat(id).getCurrentOrbitNumber();
			
			% Get updated simulation status: 0 = Stop; 1 = Good.
			[sim_status, msg] = sim.getStatus(currentOrbit);
			
			% Get time now.
			timeNow = posixtime(datetime('now')); % Posixtime [s].
			
			% Log.
			send(dq,['Sat ',num2str(id),' - Orbit ',...
			         num2str(currentOrbit),' - Duration ',...
							 num2str(timeNow - startTimeOrbit),' s - ',msg]);
			
		end % Orbital loop [while sim_status].
		
		% If orbits are broken, also break alive loop; this will change 
		% later with other conditions.
		if sim_status == 0
			alive = false;
		end
		
	end % [while alive].
	
	% Log.
	send(dq,['Sat ',num2str(id),' is dead.']);
	
end % [spmd(iv.Ns)].

% Terminate the existing parallel pool session.
delete(gcp('nocreate'));


%% Custom objects and classes used

% Set MATLAB classes to ignore.
classesToIgnore = {'Composite',...
                   'parallel.Pool',...
                   'parallel.pool.DataQueue'};

% Save current MATLAB workspace variables.
warning off parallel:lang:spmd:CompositeSave;
save(fullfile(filepath,'workspace.mat'));

% Get variables from saved workspace.
varsWorkspace = who('-file',fullfile(filepath,'workspace.mat'));
varsLength = length(varsWorkspace);

% Set empty cell array for object names.
% objectNames = {};
objectNames = cell(varsLength,1);

% Set empty cell array for class names.
% classNames = {};
classNames = cell(varsLength,1);

% Set counter for number of custom objects found.
objCounter = 0;

% Go through all variables, one by one.
for varnum = 1 : varsLength
	
	objectName = varsWorkspace{varnum};
	className = class(eval(objectName));
	
	% Check if variable is a class object.
	if isobject(eval(objectName)) && ...
	   ~any(strcmp(classesToIgnore,className))
		
		% Increment counter.
		objCounter = objCounter + 1;
		
		% Add object and class names to cell array.
		objectNames{objCounter} = objectName;
		classNames{objCounter} = className;
		
	end
	
end

if objCounter > 0
	
	% Remove empty cells from cell arrays.
	index = cellfun(@isempty,objectNames) == 0;
	objectNames = objectNames(index);
	classNames = classNames(index);
	
	% Print variable names and their respective classes.
	fprintf('\nCustom objects and respective classes:\n\n');
	objectNames = pad(objectNames,'left');
	for n = 1 : objCounter
		fprintf('%s : %s\n',objectNames{n},classNames{n});
	end
	
	% Save only class names into file.
	fid = fopen(fullfile(filepath,'lib','listCustomClasses'),'w');
	fprintf(fid,'%s\n',classNames{:});
	fclose(fid);
	
else
	fprintf('\nNo custom object classes found.\n');
end
























%% question

% initial idx and altitude
idx=120;


%iv = IvanovFormationFlight();

%  Data that will later be per satellite and therefore inside SPMD loop

orbitSection      = 2;                      % size of orbit section [deg]
orbitSectionSize  = 2;                      % size of orbit section [deg]
orbitSections     = 1:orbitSectionSize:360;
orbitCounter      = 0;
%error             = zeros(6,iv.Ns);
sst               = zeros(9,1);
sstDesired        = zeros(6,1);
sstOld            = zeros(9,1);
refPosChangeTemp  = zeros(3,1);

%  Plotting variables (will not be used for operational software)

SST_PP          = zeros(9,1); % satellite state
REFPOSCHANGE_PP = zeros(3,1); %
TIME_PP         = 0;          % time for post processing

%  Non-gravitational perturbations

% wind     = iv.WindOn * iv.Rho/2 * iv.V^2 * [-1 0 0]';
% sunlight = iv.SunOn * 2 * 4.5e-6 * [0 -1 0]'; % only for dawn/dusk orbit
% refSurf  = iv.PanelSurface * iv.Panels(3);

%  Force vector determination and angular granulaty

% alphas = 0 : iv.DeltaAngle : 360; % roll
% betas  = 0 : iv.DeltaAngle : 180; % pitch
% gammas = 0 : iv.DeltaAngle : 360; % yaw

%  Calculates pressure forces and returns 4D-Arrays of size:
%  (3, length(alphas), length(betas), length(gammas) )

% aeroPressure = aeroPressureForce(wind, iv.PanelSurface, iv.Panels(1), ...
% 	iv.Panels(2), iv.Panels(3), alphas, betas, gammas, iv.Rho, iv.V, iv.Tatmos);
% 
% solarPressure = solarPressureForce(sunlight, iv.PanelSurface, iv.Panels(1), ...
% 	iv.Panels(2), iv.Panels(3), alphas, betas, gammas);

%  Simulation object and parameters

% sim = Simulation();
% sim.MaxOrbits = 10;

%  Creates array of IvanovSatellite objects

% iv.Satellites = IvanovSatellite.empty(iv.Ns,0);
% sat = iv.Satellites; % alias for iv.Satellites
% for i = 1 : iv.Ns
% 	sat(i) = IvanovSatellite();
% end



disp('end');
