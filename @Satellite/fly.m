function fly(this, currentOrbitSection, sizeOrbitSection)
%% Initialize flight control
% ______________________________________________________________________________
%
% Details here.
% method of class Satellite
% ______________________________________________________________________________

%% rotate sunlight (note: rotation for ecliptic already done)
this.FlightControl.SolarPressure = this.FlightControl.rodriguesRotation(this.FlightControl.initialSolarPressure,[0 1 0]',currentOrbitSection/180*pi);

totalPressureDirection=this.FlightControl.WindPressure+this.FlightControl.SolarPressure;

% compute for each roll, pitch and yaw angle the solar radiation force
this.FlightControl.SolarPressureVector = this.FlightControl.getSolarPressureVector( ...
                                             this.FlightControl.SolarPressure, this.FlightControl.SurfacePanel, ...
 				                                     this.FlightControl.Panels(1), this.FlightControl.Panels(2), this.FlightControl.Panels(3), ...
 				                                     this.FlightControl.Alphas, this.FlightControl.Betas, this.FlightControl.Gammas);

% Vector of size 3 x sizeAlphas x sizeBetas x sizeGammas.
usedTotalForceVector = zeros(3,...
  size(this.FlightControl.Alphas, 2),...
  size(this.FlightControl.Betas , 2),...
  size(this.FlightControl.Gammas, 2));

masterSatellite = 0;
if masterSatellite == 0
% If no master satellite.
	for k = 1 : size(this.FlightControl.Gammas,2)
		for j = 1 : size(this.FlightControl.Betas,2)
			for i = 1 : size(this.FlightControl.Alphas,2)
				usedTotalForceVector(:,i,j,k) = this.FlightControl.WindPressureVector(:,i,j,k) + this.FlightControl.SolarPressureVector(:,i,j,k);
			end
		end
	end
else
	% To do: Implement cases with master satellite.
end


%% to-be-double-checked theory: if R is large then the control error is secondary to the minimization of the control action
R=diag([1e12 1e12 1e12]);
[P, IR, A, B] = this.FlightControl.riccatiequation(this.Orbit.MeanMotionRad, this.FlightControl.SSCoeff,R);

% determine time elapsed since last ascending equator crossing
timeSinceEqCrossing = currentOrbitSection / this.Orbit.MeanMotionDeg;

% Update desired state for all satellites.
%for i=1:this.FlightControl.NumSatellites
  this.FlightControl.updateStateDesired(timeSinceEqCrossing, this.Orbit.MeanMotionRad,i);
%end

% Set and get state error for this satellite only within array of errors of all satellites
stateError = this.FlightControl.getStateError();

%%JT: this should go through COM module% Communicate new state error of this satellite to other satellites
this.broadcastSend(stateError);
% Receive the state errors from other satellites
receivedStateErrors = this.broadcastReceive();

% Update information on state errors as received from other satellites.
this.FlightControl.updateStateErrors(receivedStateErrors);

deltaTime = sizeOrbitSection / this.Orbit.MeanMotionDeg;

this.FlightControl.StateOld = this.FlightControl.State;

oldAlphas = this.FlightControl.State(7); %% roll
oldBetas  = this.FlightControl.State(8); %% pitch
oldGammas = this.FlightControl.State(9); %% yaw

% Compute control vector for all satellites.
for i=1:this.FlightControl.NumSatellites
  this.controlVector(i,:) = (-IR * B' * P * this.FlightControl.StateErrors(:, i))';
end

errorRotation      = vrrotvec(totalPressureDirection,[1 0 0]);
angleErrorRotation = errorRotation(4)/pi*180;
axisErrorRotation  = errorRotation(1:3);

%% compute transformed control force for all satellites
for i=1:this.FlightControl.NumSatellites %% transform error for each satellite
  controlVectorTransformed(i,:)  = this.FlightControl.rodriguesRotation(this.controlVector(i,:)',axisErrorRotation',angleErrorRotation/180*pi);
end %%

%% shift and average control force for all satellites
maxControlVectorTransformed         = max(controlVectorTransformed(:,1));
for i=1:this.FlightControl.NumSatellites %% transform error for each satellite
  controlVectorTransformed(i,1)     = controlVectorTransformed(i,1) - maxControlVectorTransformed;
end
averageControlVectorTransformed     = zeros(3,1);
for i=1:this.FlightControl.NumSatellites %% compute average error
  averageControlVectorTransformed(2:3)= averageControlVectorTransformed(2:3)+controlVectorTransformed(i,2:3)'/this.FlightControl.NumSatellites;
end
for i=1:this.FlightControl.NumSatellites %% assign average error
  controlVectorTransformed(i,2:3) = controlVectorTransformed(i,2:3)-averageControlVectorTransformed(2:3)';
end

%% retransform control force for all satellites
for i=1:this.FlightControl.NumSatellites %% transform error for each satellite
  this.controlVector(i,:)              = this.FlightControl.rodriguesRotation(controlVectorTransformed(i,:)',axisErrorRotation',-angleErrorRotation/180*pi);
end 
   
%% find force vector for this satellite only    
if norm(this.controlVector(this.FlightControl.SatID,:))==0
  this.forceVector(this.FlightControl.SatID,:) = [0 0 0]'; alphaOpt=0; betaOpt=0; gammaOpt=0;
else
  [this.forceVector(this.FlightControl.SatID,:), alphaOpt, betaOpt, gammaOpt] = this.FlightControl.findBestAttitude(...
        usedTotalForceVector, this.controlVector(this.FlightControl.SatID,:),...
        this.FlightControl.Alphas, this.FlightControl.Betas, this.FlightControl.Gammas, ...
        oldAlphas, oldBetas, oldGammas);
  %% alpha is roll, beta is pitch, gamma is yaw
end
if 2*norm(this.controlVector(this.FlightControl.SatID,:))<norm(this.forceVector(this.FlightControl.SatID,:))
    this.forceVector(this.FlightControl.SatID,:)=[0 0 0]'; alphaOpt=0; betaOpt=0; gammaOpt=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vehicle translational dynamics, this part will not be used in flight software for this satellite
this.FlightControl.State(1:6) = (A * this.FlightControl.StateOld(1:6) + B * ...
                                this.forceVector(this.FlightControl.SatID,:)' / this.FlightControl.SatelliteMass) *...
                                deltaTime + this.FlightControl.StateOld(1:6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

this.FlightControl.State(7:9) = [alphaOpt betaOpt gammaOpt]'; %% roll, pitch, yaw


% Update duration of the current orbit.
this.Orbit.updateOrbitDuration();

end % Function Satellite.fly

%{
function [dSSTdt intermediateControlVector2 intermediateTime3]=myODE(intermediateTime,intermediateSST,IR,B,P,A,StateOld,this)

  persistent intermediateControlVector;
  persistent intermediateTime2;
   
  if isempty(intermediateTime)
    intermediateControlVector2=intermediateControlVector(:,2:end);
    intermediateTime3=intermediateTime2(2:end);
    dSSTdt=0;
    clear intermediateControlVector;
    clear intermediateTime2;
  else
    if isempty(intermediateControlVector)
      intermediateControlVector=zeros(3,1);
      intermediateTime2=0;
    end
    this.FlightControl.updateStateDesired(intermediateTime, this.Orbit.MeanMotionRad);
    intermediateControlVectorTemp = -IR*B'*P * (StateOld-this.FlightControl.StateDesired);
    dSSTdt=(A * StateOld + B * intermediateControlVectorTemp);
    intermediateControlVector=[intermediateControlVector intermediateControlVectorTemp];    
    intermediateTime2=[intermediateTime2 intermediateTime];
  end
  
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%
if 0%labindex==2 || labindex==3%!solve ODE with many small steps until end of timestep
  
  intermediateTime=0;
  intermediateSST=0;
  [dSSTdt, intermediateControlVector2,intermediateTime]=ode15s(@(intermediateTime,intermediateSST) myODE(intermediateTime,intermediateSST,IR,B,P,A,this.FlightControl.StateOld(1:6),this),[0 deltaTime],this.FlightControl.StateOld(1:6));
  [~, intermediateControlVector2,intermediateTime]=myODE([],[],[],[],[],[],[],[]);
  
  % compute mean of control vector
  this.controlVector(2,:)=[0 0 0];
  for i=1:3
    for j=2:size(intermediateTime,2)%% check indices in intermediateControlVector2
      this.controlVector(2,i)=this.controlVector(2,i)+intermediateControlVector2(j,i)'*(intermediateTime(j)-intermediateTime(j-1))/deltaTime;
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%}
