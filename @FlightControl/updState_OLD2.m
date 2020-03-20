function updState_OLD2(this, P, IR, A, B, deltaTime)
%% Summary here.
%_____________________________________________________________________
%
% Details here.
% method of FlightControl
%_____________________________________________________________________

this.StateOld = this.State;

oldAlphas = this.State(7);
oldBetas  = this.State(8);
oldGammas = this.State(9);

% Vector of size 3 x sizeAlphas x sizeBetas x sizeGammas.
usedTotalForceVector = zeros(3,size(this.Alphas, 2),...
	size(this.Betas , 2),size(this.Gammas, 2));

% Compute control vector.
e = this.StateErrors(:, this.SatID);
controlVector = -IR * B' * P * e;

% Sunlight rotation will be important for non dusk-dawn orbits.
rotatedSolarPressureVector = this.SolarPressureVector;

% To do:
% Rotate solar pressure vector.

masterSatellite = 0;

% If no master satellite.
if masterSatellite == 0
	for k = 1 : size(this.Gammas,2)
		for j = 1 : size(this.Betas,2)
			for i = 1 : size(this.Alphas,2)
				usedTotalForceVector(:,i,j,k) = this.WindPressureVector(:,i,j,k) + ...
                                        rotatedSolarPressureVector(:,i,j,k);
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

[forceVector, alphaOpt, betaOpt, gammaOpt] = this.findBestAttitude(...
	usedTotalForceVector, controlVector,this.Alphas, this.Betas, this.Gammas, ...
	oldAlphas, oldBetas, oldGammas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Re-check parenthesis: is that right? If so, it can be simplified.

% Update satellite state: solve ODE with backward Euler step.
this.State(1:6) = (A * this.StateOld(1:6) + B * forceVector / this.SatelliteMass) *...
                  deltaTime + this.StateOld(1:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RE-CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

this.State(7:9) = [alphaOpt betaOpt gammaOpt]';

end % Function FlightControl.updState_OLD2.
