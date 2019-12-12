%% Ivanov Formation Flight class
%___________________________________________________________________________________________________
%
%  Priority:
%  - this line 95: Are SMA and inclination equations right?
%___________________________________________________________________________________________________
%
%  To do at IvanovFormationFlightInitial.m:
%  - line 44 - check sstDesiredFunction and put IvanovFormationFlightDesired.m as function of this
%  class
%  - log.info('Original atmospheric density (rho = %1.3e) is overwritten by Ivanov\'s value (rho = 1.3e)',rho,1e-11)
%___________________________________________________________________________________________________
%
%  Done recently:
%  - orbitalproperties.m line 89: checked meaning of sma (semi-major axis)
%  - IvanovFormationFlightInitial.m line 65: checked usage of r0: it used in orbitalproperties.m for
%  calculating V
%___________________________________________________________________________________________________
%
%  References:
%  - create doc for classes: https://nl.mathworks.com/help/matlab/matlab_prog/create-help-for-classes.html
%___________________________________________________________________________________________________

classdef IvanovFormationFlight
	% IvanovFormationFlight Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		
		TotalTime = 50*90*60      % Simulation period [s], i.e approximate multiples of orbit periods (not relevant for cosmosFS)
		CompStep = 9              % Computational step size [s]
		LengthControlLoop = 90    % [s]
		TimeTemp                  % Duration and interpolation timestep for each control loop
		WakeAerodynamics = 0      % Use of wake aerodynamics
		MasterSatellite = 0       % If 0, no master; if 1, active master; if 2, passive master
		
		%
		WindOn = 1                %
		SunOn = 0                 %
		DeltaAngle = 30           % Roll, pitch, yaw angles resolution
		Ns = 4                    %
		SatelliteMass = 1         %
		Altitude = 340e3          % [m]
		ArgumentOfPerigeeAtTe0    % Not used yet
		TrueAnomalyAtTe0          % Not used yet
		EjectionVelocity = 0.5    % [m/s]
		TimeBetweenEjections = 10 % [s]
		PanelSurface = 0.01       % [m^2]
		
		Rho = 1e-11               % Constant atmospheric density [kg/m^3]
		Tatmos = 900              % This should be interpolated later [K]
		V                         % [m/s]
		RadiusOfEarth = 6371e3    % [m]
		Mu = 3.986004418e14       % Earth's standard gravitational parameter [m^3/s^2]
		MeanMotion                % [rad/s]
		R0                        % [m]
		Inclination               % [deg]
		
		SSCoeff = 1               %
		Panels = [0 0 2]          %
		SSTTemp                   %
		
		MeanAnomalyOffSet = pi/4  % For pi/2, satellite crosses on the poles; for 0, satellite crosses at equator
		NumberOfModes = 10        %
		SSParameters              %
		
		J2 = 0.00108263           % To be checked, from matlabfunctions/orbitalproperties.m
		SemiMajorAxis             % []
		
	end
	properties (Access = private)
		
		PrivatePropertyTest
		
	end
	methods
		
		function obj = IvanovFormationFlight()
			% Constructor for IvanovFormationFlight class
			
			obj.TimeTemp = 0 : obj.CompStep : obj.LengthControlLoop;
			obj.SSTTemp = zeros(9, obj.Ns, size(obj.TimeTemp,2));
			for i = 1 : obj.Ns
				obj.SSTTemp(1,i,1) = (i-1) * obj.EjectionVelocity * obj.TimeBetweenEjections; % x position
				obj.SSTTemp(4,i,1) = 0;                                                       % u velocity
				obj.SSTTemp(7,i,1) = 0;                                                       % alpha
				obj.SSTTemp(8,i,1) = 0;                                                       % beta
				obj.SSTTemp(9,i,1) = 0;                                                       % gamma
			end
			obj.SSParameters = zeros(6, obj.Ns, obj.NumberOfModes);
			obj.R0 = obj.RadiusOfEarth + obj.Altitude;
			obj.V = sqrt( obj.Mu / obj.R0 );
			obj.MeanMotion = sqrt( obj.Mu / obj.R0^3 );
			obj.SemiMajorAxis = ( 1 / obj.MeanMotion^2 * obj.Mu )^(1/3);
			obj.Inclination = acosd( -(obj.SemiMajorAxis/12352000)^(7/2) );
			
		end
		
% 		function [sstDesired] = IvanovFormationFlightDesired()
% 			
% 		end
		
	end
end