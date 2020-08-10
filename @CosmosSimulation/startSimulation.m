function startSimulation(this)
%% Initiate Cosmos simulation
% ______________________________________________________________________________
%
% Details here.
% ______________________________________________________________________________

% Create data queue for parallel pool.
dq = parallel.pool.DataQueue;

% Define function to call when new data is received on the DataQueue.
afterEach(dq, @disp);

% Create and return pool with the specified number of workers.
parpool(this.NumSatellites);

% Set the start time for the parallel pool.
timeStartPool = posixtime(datetime('now')); % Posixtime [seconds].

% Execute parallel code on workers of parallel pool.
% For better debugging, comment spmd command and its end line.
% Set id = 1 (instead of labindex).
spmd(this.NumSatellites)
  
  %! JT: most of what is done here in the parallel loop needs to go to Satellite.fly
  %! JT: process in runCosmosBeta could go to the constructor of CosmoSimulation or
  % CosmosSimulation.start. Then, maybe runCosmosBeta could be become a function of
  % CosmosSimulation. Maybe we could get rid of the aliases(they are neat from programming pov,
  % but confusing sometimes
  
	% Get unique IDs for each of the satellites, from 1 to N.
  %! JT: do we need id? isnt labindex enough?
	% for the fsw, the sat needs to find its own id in a different way, for instance from the file that read later
  id = labindex;
	
	% Create local aliases for the class objects.
	sat   = this.Satellites(id);
	orbit = this.Orbits(id);
	fc    = this.FlightControlModules(id);
  gps   = this.GPSModules(id);
  
  %!JT: we should have
  % - a class COM where we hide all details of the communication between the
  %   satellites. Communication with the screen could be kept as is.
  % - a class for communication with the screen
  % - a class for the communication with VIZ
  % Set satellite communication channel as the parpool data queue.
	commChannel = dq;
	
%!RW: transfer configuration of initial conditions to simulation, leave
%satellite initialization only with real-case-like instructions.
  % Set up some parameters, such as battery status, sat status, initial conditions.
% sat.initialize(id, commChannel, this.iniConditions(id,:));
  sat.initialize(id, commChannel, this.InitConditions(id,:));
  
    %% for sim
  % for simulation output, set initial conditions
  this.updSatStatesIni(id, fc.State);
  
  % Loop for each orbit.
	while sat.Alive % Sattelites turned on, but still doing nothing.
    
    %!R:
    % Simulate GPS data here.
    
    %% Before:
    %---------------------------------------------------------------------------
    % Update orbit counter.
    % Orbit counter holds current orbit number, not total orbits completed.
    gps.incrementOrbitCounter();
    
    % Put endOfSectionsCycle into flight control.
    % Calculate endOfSectionsCycle.
%     endOfSectionsCycle = (this.IDX - 1) / this.NumOrbitSections;
    
    % From whereInWhatOrbit().
    % This implementation can bring bugs later
%     if endOfSectionsCycle % start a new section cycle
%       gps.MeanAnomalyFromAN = 0.01;
%     else
%       gps.MeanAnomalyFromAN = 120; % when starting a simulation, the s/c might be an arbitrary mean anomaly, e.g. 120
%     end
    % ^
    % Update mean anomaly from ascending node.
    % For now, simulation updates this value.
    % Later, this value will be obtained from GPS/TLE.
    % IDX into flight control?
    %this.updateIDX(gps.MeanAnomalyFromAN);
    
    % Pause #1:
    % Wait until end of orbit sections.
    %pause( (this.OrbitSections(this.IDX) - gps.MeanAnomalyFromAN) /...
    %	orbit.MeanMotionDeg / this.AccelFactor);
    %---------------------------------------------------------------------------
    %% End
    
    %% The orbit is divided into sections of few degrees size.
    % IDX tells in which section we are in
    
    % JT: add here:
    % get meanAnomalyFromAAN from GPS:
    % meanAnomalyFromAN = sat.GPSModule.getMeanAnomalyFromAN();
    % calculate current IDX: this.updateIDX(gps.MeanAnomalyFromAN);
    % pause until the end of the orbitSection(IDX):
    % pause( (this.OrbitSections(this.IDX) - gps.MeanAnomalyFromAN) / orbit.MeanMotionDeg / this.AccelFactor);
    % JT: OrbitSections and IDX should be renamed to OrbitSectionsMeanAnomalyFromAN and orbitSectionsID
    
    timeStep = this.OrbitSectionSize / orbit.MeanMotionDeg;
    % ^^^^
    %!R:
    % Will OrbitSectionSize and MeanMotionDeg be always constant?
    
    
    %% Start orbit sections loop.
		while this.OrbitSectionNow <= this.NumOrbitSections
			
			% Determine start time of this cycle, in order to subtract the 
			% total cycle duration from the pause time (pause #2).
			timeStartSection = now();
			
			% Start flying on orbital loop.
      % in which section (in degree) are we?
%!RW: this.OrbitSectionNow holds the ID of the current orbit section, while the
%variable currentOrbitSection holds the value in degrees. This implementation is
%very confusing. Maybe change name of variables or the implemenation method to
%be more clear and avoid bugs/errors later.
      currentOrbitSection = this.OrbitSections(this.OrbitSectionNow);
      %send(DQ,'bef');
      
      % run the formation flight algorithm
%!RW: same from above notes applies here
%!RW: later, the handling of the OrbitSectionSize needs to be implemented in
%flight control module
      sat.fly(currentOrbitSection, this.OrbitSectionSize);
      %send(DQ,'after');
			
      %Before:::
			%timestep = this.OrbitSectionSize / orbit.MeanMotionDeg;
      
      % Send reference position to all non 1-satellites.
      %! this should go into a new COM module
%!RW: this implementation is a type of master-slave communication. Only 
%satellite 1 is sending the position changes to other satellites. This should be
%implemented inside of a new communication module/object.
      refPosChange = zeros(3,1);
      if id == 1
        refPosChange(1:3) = fc.State(1:3) - fc.StateOld(1:3);
        for satID = 2 : this.NumSatellites
          tag = 1000000 * satID + ...
                  10000 * this.OrbitSectionNow + ...
                    100 * orbit.OrbitCounter + ...
                      1;
          labSend(refPosChange, satID, tag);
        end
      end
      
      % Receive reference position in other satellites.
%!RW: from notes above, this is a type of master-slave communication. All
%satellites, with exception of satellite 1, are receiving the tag sent by the
%master (satellite 1). Also implement this part inside of a new communication
%module/object.
      satID = id;
      if satID ~= 1
        tag = 1000000 * satID + ...
                10000 * this.OrbitSectionNow + ...
                  100 * orbit.OrbitCounter + ...
                    1;
        refPosChange = labReceive(1, tag);
      end
      
			% Move coordinate system.
			% Should the old state be shifted as well?
			shift = -refPosChange(1:3);
      fc.shiftState(shift);
      
      % Update vector with satellite positions
      this.updSatPositions(id, refPosChange);
      % Update vector with satellite states
      this.updSatStates(id, fc.State);
      % Update time vector
      this.updTimeVector(id, timeStep);
      % add instantaneous controlVector to controlVectorTM
      sat.updControlVectorTM(id);
      % add instantaneous forceVector to forceVectorTM
      sat.updForceVectorTM();
      
%% Move to flight control
			% Increment section counter.
			this.incrementIDX();
			
			% Pause #2:
			% Pause after subtracting this section's computing time.
			pause(this.OrbitSectionSize / orbit.MeanMotionDeg /this.AccelFactor - (now() - timeStartSection));
			
		end % While orbit sections loop.
		
		% Check if orbit counter identifiers do not match.
    %JT: what is this good for?
		if (orbit.OrbitCounter ~= orbit.TimeOrbitDuration(1))
			msg = ['Orbit identifiers in orbit.OrbitCounter and ',...
				'orbit.TimeOrbitDuration do not match.'];
			error('Simulation:start:orbitIdentifierNotEqual',msg);
		else
			msg = ['Orbit ',num2str(orbit.OrbitCounter),' finished ',...
				'(',num2str(orbit.TimeOrbitDuration(2)),' s)'];
			sat.comm(msg);
      this.setIDX(1);
		end
		
    % If maximum number of orbits has been reached, turn off the satellite.
		if orbit.OrbitCounter >= this.MaxNumOrbits
% 			pause(2);
			send(dq,['[sim] Maximum number of orbits reached! ','Killing [',sat.Name,']']);
			sat.turnOff();
		end
		
	end % While alive (main orbital loop).
  
  
  % write telemetry to files, should be a function of satellite
  % also the data containers need to be properties of the Satellite, not Simulation
  writematrix(this.TimeVector(id,:)',strcat('TMTimeVector',num2str(id),'.csv'));
  writematrix(squeeze(this.SatPositions(id,:,:))',strcat('TMSatPosition',num2str(id),'.csv'));
  writematrix(squeeze(this.SatStates(id,:,:))',strcat('TMSatStates',num2str(id),'.csv'));
  writematrix(sat.controlVectorTM,strcat('TMControlVector',num2str(id),'.csv'));
  writematrix(sat.forceVectorTM,strcat('TMForceVector',num2str(id),'.csv'));
	
  % Needed for autonomous documentation generation tool.
	% Globally concatenate all output variables on lab index 1.
	% Must be the last lines of code of the parallel pool.
	satellites = gcat(sat,1,1);
	orbits = gcat(orbit,1,1);
	flightControlModules = gcat(fc,1,1);
	gpsModules = gcat(gps,1,1);
	timeVectorLengths = gcat(this.TimeVectorLengths(id),1,1);
	timeVector = gcat(this.TimeVector(id,:),1,1);
	satPositionsLengths = gcat(this.SatPositionsLengths(id),1,1);
	satPositions = gcat(this.SatPositions(id,:,:),1,1);
	satStatesLengths = gcat(this.SatStatesLengths(id),1,1);
	satStates = gcat(this.SatStates(id,:,:),1,1);
	
end % Parallel code.

% Needed for autonomous documentation generation tool.
% Get the globally concatenated values stored on lab index 1.
% Must be placed right after the end of the parallel pool.
this.Satellites = satellites{1};
this.Orbits = orbits{1};
this.FlightControlModules = flightControlModules{1};
this.GPSModules = gpsModules{1};
this.TimeVectorLengths = timeVectorLengths{1};
this.TimeVector = timeVector{1};
this.SatPositionsLengths = satPositionsLengths{1};
this.SatPositions = satPositions{1};
this.SatStatesLengths = satStatesLengths{1};
this.SatStates = satStates{1};

% Terminate the existing parallel pool session.
delete(gcp('nocreate'));

% Calculate the execution time of the parallel pool.
timeEndPool = posixtime(datetime('now')); % Posixtime [seconds].
timeDurationPool = timeEndPool - timeStartPool;
fprintf('Total simulation time: %s seconds.\n',...
num2str(timeDurationPool));

end % Function CosmosSimulation.startSimulation
