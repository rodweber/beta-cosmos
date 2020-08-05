%% Create simulation environment for satellite formation flight
% ______________________________________________________________________________
%
% Class CosmosSimulation:
%
% Detailed explanation goes here.
% ______________________________________________________________________________

classdef CosmosSimulation < handle
  
  properties (GetAccess = public, SetAccess = private)
    
    param
    iniConditions %% initial conditions
    vizScale
    
    AccelFactor % Acceleration factor for the simulation.
    FlightControlModules % Array of FlightControl objects.
    GPSModules % Array of GPS objects.
    IDX % Change this ????????????????????????????????????????????????
    MaxNumOrbits % Maximum number of orbits to run.
    NumOrbitSections % Total number of orbit sections.
    NumSatellites % Total number of satellites in the formation.
    Orbits % Array of Orbit objects.
    OrbitSectionSize % Size of each orbit section [deg].
    OrbitSections % Orbital sections for the simulation.
    Satellites % Array of Satellite objects.
    SatPositions % Satellite positions in relation to the reference.
    SatPositionsLengths % Length of the satellite positions vectors.
    SatStates % Satellites states for plotting.
    SatStatesLengths % Length of the satellite states vectors.
    SimParams % Struct holding all simulation parameters set by the user.
    %Status % Simulation status.
    TimeVector % Time vector for plotting.
    TimeVectorLengths % Length of the time vector for each satellite.
    
  end
  
  methods % Constructor.
    
    function this = CosmosSimulation(param, iniConditions)
      %% Constructor for class CosmosSimulation
      %
      % Input:
      % - Struct of parameters.
      %
      % Output:
      % - Object of class Simulation.
      
      this.param = param;
      this.iniConditions=iniConditions;
      this.vizScale = param.vizScale;
      
      this.MaxNumOrbits = param.MaxNumOrbits;
      this.AccelFactor = param.AccelFactor;
      this.IDX = param.InitIDX;
      this.OrbitSectionSize = param.OrbitSectionSize;
      this.OrbitSections = param.OrbitSectionSize:param.OrbitSectionSize:360;
      this.NumOrbitSections = length(this.OrbitSections);
      
      this.NumSatellites = param.NumSatellites;
      
      this.SatPositions = zeros(this.NumSatellites,3,1);
      this.SatPositionsLengths = ones(this.NumSatellites,1);
      
      this.SatStates = zeros(this.NumSatellites,9,1);
      this.SatStatesLengths = ones(this.NumSatellites,1);
      
      this.TimeVector = zeros(this.NumSatellites,1);
      this.TimeVectorLengths = ones(this.NumSatellites,1);
      
      % Create array with objects of class Satellite.
      this.Satellites = Satellite.empty(this.NumSatellites,0);
      
      % Get location in which to save file with FFPS for the satellites.
      ffpsFolder = param.FolderFFPS;
      [~,~,~] = mkdir(ffpsFolder); % [status,msg,msgID]
      
      % Set common part of the name for FFPS files.
      ffpsFileName = 'fc_FFP_sat';
      
      for i = 1 : this.NumSatellites
        % Create a JSON file for each satellite's formation flight parameters.
        ffpsPath = strcat(ffpsFolder,filesep,ffpsFileName,num2str(i),'.json');
        fid = fopen(ffpsPath,'w');
        fprintf(fid,'%s',jsonencode( param.FFPS(i) ) );
        fclose(fid);
        
        % Bring each satellite to life.
        % Inform the location of the FFPS file to each satellite.
        this.Satellites(i) = Satellite( ...
          param.Altitude, ...
          param.DeltaAngle, ...
          param.AutoResponse, ...
          param.AvailableGPS, ...
          param.AvailableTLE, ...
          param.NumSatellites, ...
          param.FormationMode,...
          ffpsPath);
      end
      
      % Create aliases for satellite orbits.
      this.Orbits = Orbit.empty(this.NumSatellites,0);
      for i = 1 : this.NumSatellites
        this.Orbits(i) = this.Satellites(i).Orbit;
      end
      
      % Create aliases for satellite flight control modules.
      this.FlightControlModules = FlightControl.empty(this.NumSatellites,0);
      for i = 1 : this.NumSatellites
        this.FlightControlModules(i) = this.Satellites(i).FlightControl;
      end
      
      % Create aliases for GPS modules.
      this.GPSModules = GPS.empty(this.NumSatellites,0);
      for i = 1 : this.NumSatellites
        this.GPSModules(i) = this.Satellites(i).GPSModule;
      end
      
    end % Constructor function.
    
  end % Constructor section.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Public Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  methods (Access = public)
    
    function updTimeVector(this, satID, timestep)
      % Increment length of the time vector for the current satellite.
      % At the end of the simulation, the values in TimeLength
      % should be equal for all satellites. Later, use this value to
      % prealocate memory for the TimePlot vector in order to reduce
      % computational time.
      lastPos = this.TimeVectorLengths(satID);
      nextPos = this.TimeVectorLengths(satID) + 1;
      this.TimeVector(satID, nextPos) = this.TimeVector(satID, lastPos) + timestep;
      this.TimeVectorLengths(satID) = this.TimeVectorLengths(satID) + 1;
    end
    
    function updSatStates(this, satID, satState)
      nextPos = this.SatStatesLengths(satID) + 1;
      this.SatStates(satID, 1:9, nextPos) = satState;
      this.SatStatesLengths(satID) = this.SatStatesLengths(satID) + 1;
    end
    
    function updSatStatesIni(this, satID, satState)
      this.SatStates(satID, 1:9, 1) = satState;
    end
    
    %! give better name this is the reference position change
    function updSatPositions(this, satID, newValue)
      nextPos = this.SatPositionsLengths(satID) + 1;
      this.SatPositions(satID, 1:3, nextPos) = newValue;
      this.SatPositionsLengths(satID) = this.SatPositionsLengths(satID) + 1;
    end
    
    updateIDX(this, meanAnomalyFromAN)
    startSimulation(this)
    incrementIDX(this)
    
    function setIDX(this, value)
      this.IDX = value;
    end
    
    function setSimParams(this, simParamStruct)
      this.SimParams = simParamStruct;
    end
    
  end % Public methods.
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Static)
    
    %visualizationLONLATALT(ns,ttime,sstx,ssty,sstz,pitch,yaw,roll,altitude)
    visualizationLONLATALT(ns,VIZaltitude)
    
		%plotting(angles, sst, refPosChange, time, ns, meanMotion, u, e)
    plotting(ns, meanMotionRad)
		createListCustomClasses(filepath, workspaceFileName)
		
	end % Static methods.
	
end % Class CosmosSimulation.
