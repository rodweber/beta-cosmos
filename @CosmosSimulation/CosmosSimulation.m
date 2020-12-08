%% Create simulation environment for satellite formation flight
% ______________________________________________________________________________
%
% Class CosmosSimulation:
%
% Detailed explanation goes here.
% ______________________________________________________________________________

classdef CosmosSimulation < handle

properties (GetAccess = public, SetAccess = public)

  AccelFactor % Acceleration factor for the simulation.
  AllParams % Struct with all parameters set by the user.
  FlightControlModules % Array of FlightControl objects.
  GPSModules % Array of GPS objects.
%!RW: deleted param idx: replaced it by OrbitSectionNow
% IDX % 
%!RW: change name? initial conditions for what?
  InitConditions % 
  MaxNumOrbits % Maximum number of orbits to run.
  ModeName % Formation flight mode name.
  NumOrbitSections % Total number of orbit sections.
  NumSatellites % Total number of satellites in the formation.
  Orbits % Array of Orbit objects.
  OrbitSections % Array of orbit sections with defined angle step [deg].
  OrbitSectionInitial % Initial orbit section ID for the simulation.
  OrbitSectionNow % Current orbit section ID.
  OrbitSectionSize % Size of each orbit section [deg].
  Satellites % Array of Satellite objects.
  SatPositions % Satellite positions in relation to the reference.
  SatPositionsLengths % Length of the satellite positions vectors.
%!RW: leave parameter status for later implementation.
% Status % Simulation status.
  TelemetryPath %
  VizScale % Scale for satellite distances in visualization.
  WindFactor %% switch on whether aerodynamics shall be simulated
	SolarFactor  %% switch on whether solar radiation pressure shall be simulated


end

  methods % Constructor.
    
    function this = CosmosSimulation()
      %% Constructor for class CosmosSimulation
      %
      % Input:
      % - Struct of parameters.
      %
      % Output:
      % - Object of class Simulation.
      
      % Read parameters from JSON file.
      filename = 'CosmosParameters.json';
      fid = fopen(filename,'r');
      params = jsondecode(fscanf(fid,'%s'));
      fclose(fid);
      
      % Get all parameters that are contained in the JSON file.
      modeID = params.SelectedFormationFlightMode;
      maxNumOrbits = params.SimMaxNumOrbits;
      accelFactor = params.SimAccelerationFactor;
      vizScale = params.VizScale;
      modeName = params.Modes(modeID).ModeName;
      numSatellites = params.Modes(modeID).NumSatellites;
%!RW: attitudeResolutionDeg -> former deltaAngle
      attitudeResolutionDeg = params.Modes(modeID).AttitudeResolutionDeg;
      initOrbitAltitude = params.Modes(modeID).InitOrbitAltitude;
      sizeOrbitSection = params.Modes(modeID).SizeOrbitSectionDeg;
      initOrbitSection = params.Modes(modeID).InitOrbitSectionID;
      autoResponse = params.Modes(modeID).AutoResponse;
      availableGPS = params.Modes(modeID).AvailableGPS;
      availableTLE = params.Modes(modeID).AvailableTLE;
      initConditions = params.Modes(modeID).InitialConditions;
      ffpsFolderName = params.Modes(modeID).FFPSFolder;
      ffpsValues = params.Modes(modeID).FFPSValues;
      
      % Set simulation parameters.
      this.AccelFactor = accelFactor;
      this.AllParams = params;
      this.MaxNumOrbits = maxNumOrbits;
      this.ModeName = modeName;
      this.NumSatellites = numSatellites;
      this.VizScale = vizScale;
      
%!RW: Parameters still in simulation, later transfer to FlightControl.
      this.InitConditions = initConditions;
      this.OrbitSections = sizeOrbitSection:sizeOrbitSection:360;
      this.OrbitSectionInitial = initOrbitSection;
      this.OrbitSectionNow = initOrbitSection;
      this.OrbitSectionSize = sizeOrbitSection;
      this.NumOrbitSections = length(this.OrbitSections);
            
      % Create array for objects of class Satellite.
      this.Satellites = Satellite.empty(numSatellites,0);
      
      % Create arrays for aliases of Satellite children objects.
      % Aliases allow calling children objects with shorter commands.
      this.Orbits = Orbit.empty(numSatellites,0);
      this.FlightControlModules = FlightControl.empty(numSatellites,0);
      this.GPSModules = GPS.empty(numSatellites,0);
      
      % Set location in which to save files with FFPS for the satellites.
      ffpsFolderPath = strcat('storage',filesep,ffpsFolderName);
      [~, ~, ~] = mkdir(ffpsFolderPath); % [status, msg, msgID]
      ffpsFilePrefix = 'fc_FFP_sat';
      ffpsPathPrefix = strcat(ffpsFolderPath,filesep,ffpsFilePrefix);
      
      % Set location for telemetry files.
      tmFolderPath = 'telemetry';
      [~, ~, ~] = mkdir(tmFolderPath); % [status, msg, msgID]
      this.TelemetryPath = tmFolderPath;
      
      for i = 1 : numSatellites
        % Create a JSON file for each satellite's formation flight parameters.
        ffpsFullPath = strcat(ffpsPathPrefix,num2str(i),'.json');
        fid = fopen(ffpsFullPath,'w');
        fprintf(fid,'%s',jsonencode( ffpsValues(i) ) );
        fclose(fid);
        
        % Bring each satellite to life.
        % Inform the location of the FFPS file to each satellite.
        this.Satellites(i) = Satellite( ...
          initOrbitAltitude, ...
          attitudeResolutionDeg, ...
          autoResponse, ...
          availableGPS, ...
          availableTLE, ...
          numSatellites, ...
          modeID, ...
          ffpsFullPath, ...
          tmFolderPath);
        
        % Create simulation aliases for satellite children objects.
        this.Orbits(i) = this.Satellites(i).Orbit;
        this.FlightControlModules(i) = this.Satellites(i).FlightControl;
        this.GPSModules(i) = this.Satellites(i).GPSModule;
      end
      
      this.WindFactor = 1;      %% shall aerodynamics be simulated?
      this.SolarFactor = 1;     %% shall solar radiation pressure be simulated?
      
    end % Constructor function.
    
  end % Constructor section.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Public Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods (Access = public)
    
    updateIDX(this, meanAnomalyFromAN)
    startSimulation(this)
    incrementIDX(this)
    
%!RW: later rename function to setOrbitSection?
    function setIDX(this, value)
%       this.IDX = value;
      this.OrbitSectionNow = value;
    end
    
    %!RW: for reference, old function:
    % plotting(angles, sst, refPosChange, time, ns, meanMotion, u, e)
    plotting(this, ns, meanMotionRad)
    
  end % Public methods.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)
  
  [time,lat,lon,rad]=  keplerPropagation(cosmosTime,keplerStepSize,inclination,RAAN,v0,altitude,radiusOfEarth)
  
  %!RW: for reference, old function:
  % visualizationLONLATALT(ns,ttime,sstx,ssty,sstz,pitch,yaw,roll,altitude)
  % visualizationLONLATALT(this, vizScale, ns, altitude)
  ECEFprocessing(this,VIZscale,ns,altitude,radiusOfEarth)

  GNSSRprocessing(this,ns,radiusOfEarth)

  createListCustomClasses(filepath, workspaceFileName)
  
end % Static methods.

end % Class CosmosSimulation.
