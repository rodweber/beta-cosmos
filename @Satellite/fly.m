function fly(this, currentOrbitSection, sizeOrbitSection)
%% Initialize flight control
% ______________________________________________________________________________
%
% Details here.
% method of class Satellite
% ______________________________________________________________________________

% Get updated orbital parameters from GPS/TLE.
orbitFromGPS = this.GPSModule.getOrbitCounter();
meanAnomalyFromAN = this.GPSModule.getMeanAnomalyFromAN();

% Update orbital parameters.
this.Orbit.updateOrbitalParams(orbitFromGPS, meanAnomalyFromAN);

% Update wind pressure.
this.FlightControl.updWindPressures(this.Orbit.Rho, this.Orbit.V, this.Orbit.TempAtmos);

% Settings for control algorithm, is this necessary in every loop?
[P, IR, A, B] = this.FlightControl.riccatiequation(this.Orbit.MeanMotionRad, this.FlightControl.SSCoeff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From main cosmosFS:
% Set attitude computed in last iteration.
% this.FlightControl.setAttitude();

% Compute attitude for next section.
% Determine desired trajectory.


% determine time elapsed since last ascending equator crossing
time = currentOrbitSection / this.Orbit.MeanMotionDeg;

% Update desired state.
this.FlightControl.updateStateDesired(time, this.Orbit.MeanMotionRad);

% apply noise
% maxRandPos=0;
% maxRandVel=0.;
% this.FlightControl.State(1)=this.FlightControl.State(1)+maxRandPos*(rand-0.5);
% this.FlightControl.State(2)=this.FlightControl.State(2)+maxRandPos*(rand-0.5);
% this.FlightControl.State(3)=this.FlightControl.State(3)+maxRandPos*(rand-0.5);
% this.FlightControl.State(4)=this.FlightControl.State(4)+maxRandVel*(rand-0.5);
% this.FlightControl.State(5)=this.FlightControl.State(5)+maxRandVel*(rand-0.5);
% this.FlightControl.State(6)=this.FlightControl.State(6)+maxRandVel*(rand-0.5);

% Set and get state error for this satellite in array of errors of all satellites
stateError = this.FlightControl.getStateError();

%%JT: this should go through COM module
% Communicate new state error of this satellite to other satellites
this.broadcastSend(stateError);
% Receive the state errors from other satellites
receivedStateErrors = this.broadcastReceive();

% Update information on state errors from other satellites.
this.FlightControl.updateStateErrors(receivedStateErrors);



%{
NOT USEFUL ANYMORE SINCE WORKING WITH FORCEVECTOR AVERAGES/SHIFTS
% Calculate average of the state errors.
avg = this.FlightControl.getStateErrorAverage();

% Communicate the calculated average error to the other satellites.
this.broadcastSend(avg);

% Receive the average state errors from other satellites.
receivedAverageStateErrors = this.broadcastReceive();

% compare different average state vectors as a fault detection methods

% JT: I don't know what this is good for
% Update information on the average state errors.
%this.FlightControl.updateStateErrorsAvg(receivedAverageStateErrors);

% Shift error in main force direction
%JT: remove:this.FlightControl.shiftError();
%JT: add
this.FlightControl.StateErrors(1,this.FlightControl.SatID) = this.FlightControl.StateErrors(1,this.FlightControl.SatID) - max(this.FlightControl.StateErrors(1,:));

% average error

% From main cosmosFS:
% Check whether all errors have the same value in all satellites.
%}

% Compute attitude and update satellite state.
deltaTime = sizeOrbitSection / this.Orbit.MeanMotionDeg;
%deltaTime
%JT removed:this.FlightControl.updState(P, IR, A, B, deltaTime);

%JT added:----
%------------------------------------------------
this.FlightControl.StateOld = this.FlightControl.State;

oldAlphas = this.FlightControl.State(7);
oldBetas  = this.FlightControl.State(8);
oldGammas = this.FlightControl.State(9);

% Vector of size 3 x sizeAlphas x sizeBetas x sizeGammas.
usedTotalForceVector = zeros(3,...
  size(this.FlightControl.Alphas, 2),...
  size(this.FlightControl.Betas , 2),...
  size(this.FlightControl.Gammas, 2));

% Compute control vector.
for i=1:this.FlightControl.NumSatellites
  this.controlVector(:,i) = -IR * B' * P * this.FlightControl.StateErrors(:, i);
end

%%

    %% shift and average control force
    controlVectorMin         = min(this.controlVector(1,:));
    for i=1:this.FlightControl.NumSatellites %% transform error for each satellite
      this.controlVector(1,i)     = this.controlVector(1,i) - controlVectorMin;
    end
    averagecontrolVector     = zeros(3,1);
    for i=1:this.FlightControl.NumSatellites %% compute average error
      averagecontrolVector(2:3)= averagecontrolVector(2:3)+this.controlVector(2:3,i)/this.FlightControl.NumSatellites;      
    end
    for i=1:this.FlightControl.NumSatellites %% assign average error
      this.controlVector(2:3,i) = this.controlVector(2:3,i)-averagecontrolVector(2:3);
    end
%%





% Sunlight rotation will be important for non dusk-dawn orbits.
rotatedSolarPressureVector = this.FlightControl.SolarPressureVector;

% To do:
% Rotate solar pressure vector.

masterSatellite = 0;

% If no master satellite.
if masterSatellite == 0
	for k = 1 : size(this.FlightControl.Gammas,2)
		for j = 1 : size(this.FlightControl.Betas,2)
			for i = 1 : size(this.FlightControl.Alphas,2)
				usedTotalForceVector(:,i,j,k) = this.FlightControl.WindPressureVector(:,i,j,k) + rotatedSolarPressureVector(:,i,j,k);
			end
		end
	end
else
	% To do:
	% Implement cases with master satellite.
	% Will the formation mode with master satellite impact other 
	% functions as well? Check it before implementing.
	
	% Wind and sunlight for dusk-dawn orbit.
	% wind     = norm(this.WindPressure)  * [-1 0 0]';
	% sunlight = norm(this.SolarPressure) * [0 -1 0]';
	
	% solarForceVectorOfMaster = [0 0 0]';
	% maxSolarForce            = [0 0 0]';
	% forceOnMaster = 0.5 * 2.8 * (wind + sunlight) * this.SurfaceRef;
	
	% usedTotalForceVector = ...
end

if norm(this.controlVector(:,this.FlightControl.SatID))==0
  this.forceVector = [0 0 0]'; alphaOpt=0; betaOpt=0; gammaOpt=0;
else
  [this.forceVector, alphaOpt, betaOpt, gammaOpt] = this.FlightControl.findBestAttitude(...
    usedTotalForceVector, this.controlVector(:,this.FlightControl.SatID),...
    this.FlightControl.Alphas, this.FlightControl.Betas, this.FlightControl.Gammas, ...
    oldAlphas, oldBetas, oldGammas);
end
if 2*norm(this.controlVector(:,this.FlightControl.SatID))<norm(this.forceVector)
    this.forceVector=[0 0 0]'; alphaOpt=0; betaOpt=0; gammaOpt=0;
end

method2SolveODE=1;
switch method2SolveODE %% update satellite state through solving ODE with
  case 1 %  with backward Euler step.
      this.FlightControl.State(1:6) = (A * this.FlightControl.StateOld(1:6) + B * this.forceVector / this.FlightControl.SatelliteMass) *...
                                deltaTime + this.FlightControl.StateOld(1:6);
  case 2
      ;
  otherwise
      ;
end

%% Before:
% Update satellite state: solve ODE with backward Euler step.
% this.FlightControl.State(1:6) = (A * this.FlightControl.StateOld(1:6) + B * this.forceVector / this.FlightControl.SatelliteMass) *...
%                                 deltaTime + this.FlightControl.StateOld(1:6);

this.FlightControl.State(7:9) = [alphaOpt betaOpt gammaOpt]';

%------------------------------------------------
%JT: end----






% Update duration of the current orbit.
this.Orbit.updateOrbitDuration();


end % Function Satellite.fly
