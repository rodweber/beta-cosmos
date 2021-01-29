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
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------  
spmd(this.NumSatellites)
  % For DEBUGGING mode:
  % 1. Comment the spmd command above and its end line.
  % 2. Set variable uid = 1, instead of default uid = labindex.
  
  % Set a Unique Identifier (UID) for each of the satellites, from 1 to N.
  % Create UIDs as aliases to labindex for future compatibility with C++ and to
  % increase code readability and maintainability with shorter code lines.
  uid = labindex;
  
%%TODO: for the fsw, the sat needs to find its own id in a different way, for
%%instance from the file that read later.
  
  %! JT: most of what is done here in the parallel loop needs to go to Satellite.fly
  %! JT: process in runCosmosBeta could go to the constructor of CosmoSimulation or
  % CosmosSimulation.start. Then, maybe runCosmosBeta could be become a function of
  % CosmosSimulation.
  
  % Create local aliases for the class objects.
  % Aliases are used to eliminate calls of deep functions, preventing bugs.
  sat   = this.Satellites(uid);
  orbit = this.Orbits(uid);
  fc    = this.FlightControlModules(uid);
  gps   = this.GPSModules(uid);
  
  %!JT: we should have
  % - a class COM where we hide all details of the communication between the
  %   satellites. Communication with the screen could be kept as is.
  % - a class for communication with the screen
  % - a class for the communication with VIZ
  
%!RW:TODO:
% 
  
  % Initialize satellite.
  % Set between-satellite communication channel as the parpool data queue.
  % Satellite automatically reads initial conditions from file in storage.
  sat.initialize(uid, dq, this.InitConditions(uid,:));
  %^^^^^^
%!RW: transfer configuration of initial conditions to simulation, leave
%satellite initialization only with real-case-like instructions.
%!RW: TODO:
% Simulation must read configuration files with initial conditions for satellite
% config, formation flight config, and simulation config. Satellite and
% formation flight configurations are saved in satellite storage. Simulation
% config is loaded into parameters of Simulation object.
  
  lastTime=0; %% lastTime is used to wrap around the time vector for each new orbit
  
  % define nominal wind magnitude and direction
  sat.FlightControl.WindPressure = this.WindFactor * sat.Orbit.Rho/2 * sat.Orbit.V^2 * [-1 0 0]';
  % compute for each roll, pitch and yaw angle the aerodynamic force
  sat.FlightControl.WindPressureVector = FlightControl.getWindPressureVector(...
                                              sat.FlightControl.WindPressure, sat.FlightControl.SurfacePanel, ...
	                                            sat.FlightControl.Panels(1), sat.FlightControl.Panels(2), ...
                                              sat.FlightControl.Panels(3), sat.FlightControl.Alphas, sat.FlightControl.Betas,...
                                              sat.FlightControl.Gammas, sat.Orbit.Rho, sat.Orbit.V, sat.Orbit.TempAtmos);


  % define nominal solar radiation pressure magnitude and direction
  sat.FlightControl.initialSolarPressure = this.SolarFactor * 2 * 4.5e-6 * [0 -1 0]';

  %% incline by ecliptic
  sat.FlightControl.initialSolarPressure = sat.FlightControl.rodriguesRotation(sat.FlightControl.initialSolarPressure,[0 0 1]',-23.4/180*pi);
  
  %% for sim
  % for simulation output, set initial conditions
  sat.updSatStatesIni(uid, fc.State);
  
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------   
  % Loop for each orbit.
	while sat.Alive % Satellites turned on, but still doing nothing.
        
    % Update orbit counter.
    % Orbit counter holds current orbit number, not total orbits completed.
    gps.incrementOrbitCounter();
    
    % Get updated orbital parameters from GPS/TLE.
    orbitFromGPS = sat.GPSModule.getOrbitCounter();
    meanAnomalyFromAN = sat.GPSModule.getMeanAnomalyFromAN();
    
    % Update orbital parameters.
    sat.Orbit.updateOrbitalParams(orbitFromGPS, meanAnomalyFromAN);
            
    %% ISL COM: ECEF coordinates          
    
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
    
%-------------------------------------------------------------------------------  
    %% Start orbit loop
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
      % run the formation flight algorithm
      %%% THIS SHOULD GO TO SATELLITE OR SATELLITE.FLIGHTCONTROL
      sat.fly(currentOrbitSection, this.OrbitSectionSize);
      
      %%%%%%%%THIS SHOULD GO TO AN ISL COM MODULE
      %% currently, the reference position change is communicated to the other sats
      %% this needs to be rethought towards exchange of Hill coordinates/distances
      %% sending:
%!RW: this implementation is a type of master-slave communication. Only 
%satellite 1 is sending the position changes to other satellites. This should be
%implemented inside of a new communication module/object.
      refPosChange = zeros(3,1);
      if uid == 1
        refPosChange(1:3) = fc.State(1:3) - fc.StateOld(1:3);
        for i = 2 : this.NumSatellites
          tag = 1000000 * i + ...
                  10000 * this.OrbitSectionNow + ...
                    100 * orbit.OrbitCounter + ...
                      1;
          labSend(refPosChange, i, tag);
        end
      end
      %% receiving:
%!RW: from notes above, this is a type of master-slave communication. All
%satellites, with exception of satellite 1, are receiving the tag sent by the
%master (satellite 1). Also implement this part inside of a new communication
%module/object.
      if uid ~= 1
        tag = 1000000 * uid + ...
                10000 * this.OrbitSectionNow + ...
                  100 * orbit.OrbitCounter + ...
                    1;
        refPosChange = labReceive(1, tag);
      end

			% Move coordinate system.
      %%% Should the old state be shifted as well?
      fc.shiftState(-refPosChange(1:3));            

      % Update (orbit) TM vector with satellite positions
      sat.updSatPositionsTM(uid, refPosChange);
      % Update (orbit) TM vector with satellite states
      sat.updSatStatesTM(uid, fc.State);
      % Update (orbit) TM time vector with time vector, save last time for swrapping around
      sat.updTimeVectorTM(uid, timeStep,lastTime);
      lastTime=sat.TimeVectorTM(uid, end);
      % Update (orbit) TM controlVector with controlVector
      sat.updControlVectorTM();
      % Update (orbit) TM forceVector with forceVector
      sat.updForceVectorTM();

%% Move to flight control
			% Increment section counter.
			this.incrementIDX();
			
			%% Pause processing until satellite is in next section
      %% to this end: compute time of next section and subtract this section's computing time.
			pause(this.OrbitSectionSize / orbit.MeanMotionDeg /this.AccelFactor - (now() - timeStartSection));
			
		end % While orbit loop
%-------------------------------------------------------------------------------  
		
		% Check if orbit counter identifiers do not match.
		if (orbit.OrbitCounter ~= orbit.TimeOrbitDuration(1))
			msg = ['Orbit identifiers in orbit.OrbitCounter and ',...
				'orbit.TimeOrbitDuration do not match.'];
			error('Simulation:start:orbitIdentifierNotEqual',msg);
    else
      %% write mission TM, i.e. append orbit TM to mission TM file, then reset orbit TM variables
      sat.writeAndResetMissionTM(uid)      
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
		
	end % While alive (main loop)
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------  
  % Needed for autonomous documentation generation tool.
	% Globally concatenate all output variables on lab index 1.
	% Must be the last lines of code of the parallel pool.
  satellites = gcat(sat,1,1);
  orbits = gcat(orbit,1,1);
  flightControlModules = gcat(fc,1,1);
  gpsModules = gcat(gps,1,1);
end % Parallel code.
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------  
%-------------------------------------------------------------------------------  

% Needed for autonomous documentation generation tool.
% Get the globally concatenated values stored on lab index 1.
% Must be placed right after the end of the parallel pool.
this.Satellites = satellites{1};
this.Orbits = orbits{1};
this.FlightControlModules = flightControlModules{1};
this.GPSModules = gpsModules{1};

% Terminate the existing parallel pool session.
delete(gcp('nocreate'));

% Calculate the execution time of the parallel pool.
timeEndPool = posixtime(datetime('now')); % Posixtime [seconds].
timeDurationPool = timeEndPool - timeStartPool;
fprintf('Total simulation time: %s seconds.\n',...
num2str(timeDurationPool));

end % Function CosmosSimulation.startSimulation

