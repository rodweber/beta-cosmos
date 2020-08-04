%% Instantiate flight control objects.
%_____________________________________________________________________
%
% Class FlightControl:
%
% Detailed explanation goes here.
%_____________________________________________________________________

classdef FlightControl < handle
	
	properties (GetAccess = public, SetAccess = public)
		
        DeltaAngle % Angular granularity for force vector determination.
        FFPS % Formation flight parameters.
        FFPSFilePath % Location of the file with the FFPS.
		FormationMode % Mode for the satellites formation flight.
		NumSatellites % Total number of satellites in the formation.
		Panels % Satellite panels.
		SatelliteMass % Mass of the satellite [kg].
		SatID % Unique identification number of the satellite.
		SolarFactor % Multiplication factor for solar pressure.
		SolarPressure % Pressure from sunlight.
		SolarPressureVector % Solar pressure for all satellite attitudes.
		SSCoeff % Schweighart-Sedgwick coefficient.
		State % Satellite state.
		StateDesired % Desired satellite state.
		StateErrors % State errors of all satellites in the formation.
		StateErrorsAvg % Average of the state errors of all satellites.
		StateOld % Old satellite state.
		SurfacePanel % Surface area of each of the satellite panels [m^2].
		SurfaceRef % Reference surface area.
		WindFactor % Multiplication factor for wind pressure.
		WindPressure % Pressure from wind.
		WindPressureVector % Wind pressure for all satellite attitudes.
        
  end
	
	properties (GetAccess = public, SetAccess = private)
		
		Alphas % Roll.
		Betas % Pitch.
		Gammas % Yaw.
		
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Constructor Method %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods % Constructor.
		
		function this = FlightControl(ns, mode, deltaAngle, ffpsFilePath)
%% Constructor for class FlightControl.
%_____________________________________________________________________
%
% Input:
% - Total number of satellites in the formation.
% - Mode for the satellites formation flight.
%
% Output:
% - Object of class FlightControl.
%_____________________________________________________________________
			
            this.FFPSFilePath = ffpsFilePath;
            this.NumSatellites  = ns;
			this.FormationMode  = mode;
			this.State          =[-950 0 0 0 0 0 0 0 0]';% zeros(9, 1);
			this.StateDesired   = zeros(6, 1);
			this.StateErrors    = zeros(6, ns);
			this.StateErrorsAvg = zeros(6, ns);
			
			this.SSCoeff = 1;
			
			% Force vector determination and angular granularity.
            this.DeltaAngle = deltaAngle; % Angular granularity.
			this.Alphas = 0:deltaAngle:360; % Roll.
			this.Betas  = 0:deltaAngle:180; % Pitch.
			this.Gammas = 0:deltaAngle:360; % Yaw.
			
			%this.SatelliteMass = 1; % Kilogram(s).
			this.SatelliteMass = 1.5; % Kilogram(s).
			%this.Panels = [0 0 2];
			%this.Panels = [0 0 3.5];
			this.Panels = [0 0 4];
			this.SurfacePanel = 0.01; % Squared meters.
			this.SurfaceRef = this.SurfacePanel * this.Panels(3);
			
			this.WindFactor = 1; 
			this.SolarFactor = 0;
			
			% From main cosmosFS:
			% Constant value for solar pressure.
			this.SolarPressure = this.SolarFactor * 2 * 4.5e-6 * [0 -1 0]';
			
			% Compute the experienced force caused by solar pressure for all
			% satellite attitudes (roll, pitch, yaw).
			this.SolarPressureVector = this.getSolarPressureVector( ...
				this.SolarPressure, this.SurfacePanel, ...
				this.Panels(1), this.Panels(2), this.Panels(3), ...
				this.Alphas, this.Betas, this.Gammas);
			
		end
		
	end % Constructor.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Public Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Access = public)
		
		updWindPressures(this, rho, v, tempAtmos)
		
		updateStateErrorsAvg(this, receivedAverageStateErrors)
		avg = getStateErrorAverage(this)
		
		updateStateErrors(this, receivedStateErrors)
		stateErr = getStateError(this)
		
		updateStateDesired(this, time, meanMotion)
		
		% Testing git/uml with method updState
		updStateold2(this, P, IR, A, B, deltaTime)
		
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
			alphas, betas, gammas, rho, v, Tatmos)
		
		solarpressureforcevector = getSolarPressureVector( ...
			sunlight, panelSurface, noxpanels, noypanels, nozpanels, ...
			alphas, betas, gammas)
		
		vRot = rodriguesRotation(v, k, theta)
		
		[P, IR, A, B] = riccatiequation(meanMotion, SSCoeff)
		
		[forceVector, alphaOpt, betaOpt, gammaOpt] = findBestAttitude( ...
			totalForceVector, controlVector, alphas, betas, gammas, ...
			oldAlphaOpt, oldBetaOpt, oldGammaOpt)
		
		[CD, CL] = aeroDragLiftSentman(theta, Tatmos, v, rho)
		
	end % Static methods.
	
end % Class FlightControl.
