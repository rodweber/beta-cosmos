%% Instantiate flight control objects
% ______________________________________________________________________________
%
% Class FlightControl:
%
% Detailed explanation goes here.
% ______________________________________________________________________________

classdef FlightControl < handle

properties (GetAccess = public, SetAccess = public)
  
  DeltaAngle % Angular granularity for force vector determination.
  FFPS % Formation flight parameters.
  FFPSFilePath % Location of the file with the FFPS.
  
  TelemetryPath %
  
  FormationMode % Mode for the satellites formation flight.
  NumSatellites % Total number of satellites in the formation.
  Panels % Satellite panels.
  SatelliteMass % Mass of the satellite [kg].
  SatID % Unique identification number of the satellite.
  SolarPressure % Pressure from sunlight. Will be rotated over the orbit
  initialSolarPressure % initial Pressure from sunlight.
  SolarPressureVector % Solar pressure for all satellite attitudes.
  SSCoeff % Schweighart-Sedgwick coefficient.
  State % Satellite state.
  StateDesired % Desired satellite state.
  StateErrors % State errors of all satellites in the formation.
  StateErrorsAvg % Average of the state errors of all satellites.
  StateOld % Old satellite state.
  SurfacePanel % Surface area of each of the satellite panels [m^2].
  SurfaceRef % Reference surface area.
  WindPressure % Pressure from wind.
  WindPressureVector % Wind pressure for all satellite attitudes.
  
end

properties (GetAccess = public, SetAccess = private)
  
  rollAngles % Roll.
  pitchAngles % Pitch.
  yawAngles % Yaw.
  
end
	
	methods % Constructor.
		
    function this = FlightControl(ns, mode, deltaAngle, ffpsFilePath, tmFolderPath)
      %% Constructor for class FlightControl
      %
      % Input:
      % - Total number of satellites in the formation.
      % - Mode for the satellites formation flight.
      %
      % Output:
      % - Object of class FlightControl.
      
      % Read formation flight parameters from JSON file.
      %%% Already done in Satellite.initialize
      %fid = fopen(ffpsFilePath,'r');
      %this.FFPS = jsondecode(fscanf(fid,'%s'));
      %fclose(fid);
      
      % Save path to formation flight parameters.
			this.FFPSFilePath = ffpsFilePath;
      
      % Save path to telemetry files.
      this.TelemetryPath = tmFolderPath;
      
      this.NumSatellites = ns;
			this.FormationMode = mode;
			this.StateDesired   = zeros(6, 1);
			this.StateErrors    = zeros(6, ns);
			this.StateErrorsAvg = zeros(6, ns);
			
			this.SSCoeff = 1;
			
			% Force vector determination and angular granularity.
      this.DeltaAngle = deltaAngle; % Angular granularity.
			this.rollAngles = 0:deltaAngle:360; % Roll.
			this.pitchAngles  = 0:deltaAngle:180; % Pitch.
			this.yawAngles = 0:deltaAngle:360; % Yaw.
			
			%this.SatelliteMass = 1; % Kilogram(s).
			this.SatelliteMass = 2; % Kilogram(s).
			%this.Panels = [0 0 2];
			this.Panels = [0 0 6];
			this.SurfacePanel = 0.01; % Squared meters.
			this.SurfaceRef = this.SurfacePanel * this.Panels(3);
						
		end
		
	end % Constructor.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Public Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Access = public)
				
		updateStateErrorsAvg(this, receivedAverageStateErrors)
		avg = getStateErrorAverage(this)
		
		updateStateErrors(this, receivedStateErrors)
		stateErr = getStateError(this)
		
		updateStateDesired(this, time, meanMotion)
		
		function shiftState(this, shift)
			this.State(1:3) = this.State(1:3) + shift;
		end
		
%JT: remove
%		function shiftError(this)
%			this.StateErrors(1,this.SatID) = this.StateErrors(1,this.SatID) - max(this.StateErrors(1,:));
%		end
		
		function updateSatelliteID(this, satID)
			this.SatID = satID;
		end
		
		function setFormationMode(this, formationMode)
			if ~isnumeric(formationMode)
				error('Formation flight mode must be numeric.');
			end
			this.FormationMode = formationMode;
		end
		
		function mode = getFormationMode(this)
			mode = this.FormationMode;
		end
		
	end % Public methods.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Static)
		
		aerototalforcevector = getWindPressureVector( ...
			wind, panelSurface, noxpanels, noypanels, nozpanels, ...
			rollAngles, pitchAngles, yawAngles, rho, v, Tatmos)
		
		solarpressureforcevector = getSolarPressureVector( ...
			sunlight, panelSurface, noxpanels, noypanels, nozpanels, ...
			rollAngles, pitchAngles, yawAngles)
		
		vRot = rodriguesRotation(v, k, theta)
		
		[P, IR, A, B] = riccatiequation(meanMotion, SSCoeff,R)
		
		[forceVector, rollAngleOpt, pitchAngleOpt, yawAngleOpt] = findBestAttitude( ...
			totalForceVector, controlVector, rollAngles, pitchAngles, yawAngles, ...
			oldRollAngleOpt, oldPitchAngleOpt, oldYawAngleOpt)
		
		[CD, CL] = aeroDragLiftSentman(theta, Tatmos, v, rho)
		sunforceVector=solarDragLift(solarPressure, sunlight,normal, theta,...
      panelSurface,noPanels,yawAngleSunSpecular,yawAngleSunDiffusive)
    
  end % Static methods
  
end % Class FlightControl
