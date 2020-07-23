%% Instantiate satellite objects.
%_____________________________________________________________________
%
% Class Satellite:
%
% Detailed explanation goes here.
%_____________________________________________________________________

classdef Satellite < handle
	
	properties (GetAccess = public, SetAccess = public)		
		Alive % If satellite is alive or not [true/false].
		AutoResponse % If satellite should send responses [true/false].
		CommChannel % Communication channel.
		FlightControl % Object of class FlightControl.
		GPSModule % Object of class GPS.
		Name % Unique name for identification of the satellite.
		Orbit % Object of class Orbit.
	end
	
	properties (GetAccess = public, SetAccess = private)
		%% n/a
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Constructor Method %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods % Constructor.
		
		function this = Satellite(altitude, deltaAngle, autoResponse, ...
			gpsAvailability, tleAvailability, numSats, mode)
      %% Constructor for class Satellite.
      %_____________________________________________________________________
      %
      % Input:
      % - Altitude [meters].
      % - Roll, pitch, yaw angles resolution [deg].
      % - Auto response activation [boolean].
      % - GPS availability [boolean].
      % - TLE availability [boolean].
      % - Total number of satellites in the formation.
      % - Mode for the satellites formation flight.
      %
      % Output:
      % - Object of class Satellite.
      %_____________________________________________________________________			
			this.AutoResponse = autoResponse;
			this.FlightControl = FlightControl(numSats, mode, deltaAngle);
			this.Orbit = Orbit(altitude, gpsAvailability, tleAvailability);
			this.GPSModule = GPS();	
		end
		
	end % Constructor.
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Public Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	methods (Access = public)	
		initialize(this, id, commChannel,iniConditinos)
		comm(this, msg)
		fly(this, currentOrbitSection, sizeOrbitSection)
		broadcastSend(this, msg)
		msg = broadcastReceive(this)
	end % Public methods.
	
end % Class Satellite.
