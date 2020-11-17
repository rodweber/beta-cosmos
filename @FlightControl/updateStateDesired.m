function updateStateDesired(this, time, meanMotionRad,satID)
%% Update desired satellite state
% ______________________________________________________________________________
%
% Desired solution for Ivanov's case.
%
% Input:
% - Mean motion [rad/s].
% ______________________________________________________________________________


% Does this make sense?
% Isn't the size of 'time' always 1?
% sstDesired = zeros(9,size(time,2));
sstDesired = zeros(6,1);

SSC=1;  %% SSCoefficient [-]

A=2*SSC/(2-SSC^2); %% if SSC==1 then A=2
B=(2-5*SSC^2)/2/SSC;
D=sqrt(3*SSC^2-2); %% if SSC==1 then D=1

%ffp1 for x and z sin/cos
%ffp2 x-direction constant [m]
%ffp3 for y  sin/cos [m]
%ffp4 for y  sin/cos [m]
%ffp5 for x and z  sin/cos [m]
%ffp6 for x linear
%ffp7 mean anomaly offset x,z [which dimension?]
%ffp8 mean anomaly offset y   [which dimension?]

sstDesired(1,:)= this.FFPS.ffp2 + B*meanMotionRad*this.FFPS.ffp6*time + this.FFPS.ffp1* cos( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 )                               - this.FFPS.ffp5 *sqrt(2*SSC*A)*   sin( sqrt(2*SSC/A)*meanMotionRad* time  + this.FFPS.ffp7 );
sstDesired(2,:)=                                                        this.FFPS.ffp3* cos( D*meanMotionRad*             time + this.FFPS.ffp8 )                               + this.FFPS.ffp4 /D/meanMotionRad* sin( D*meanMotionRad*             time  + this.FFPS.ffp8 );
sstDesired(3,:)= this.FFPS.ffp6                                       + this.FFPS.ffp5* cos( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 )                               + this.FFPS.ffp1 /sqrt(2*SSC*A)*   sin( sqrt(2*SSC/A)*meanMotionRad* time  + this.FFPS.ffp7 );

sstDesired(4,:)=                  B*meanMotionRad*this.FFPS.ffp6      - this.FFPS.ffp1* sin( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 ) *sqrt(2*SSC/A)*meanMotionRad  - this.FFPS.ffp5 *sqrt( 2*SSC*A)*  cos( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 ) *sqrt(2*SSC/A)*meanMotionRad;
sstDesired(5,:)=                                                      - this.FFPS.ffp3* sin( D*meanMotionRad*             time + this.FFPS.ffp8 ) *D*meanMotionRad              + this.FFPS.ffp4 /D/meanMotionRad* cos( D*meanMotionRad*             time + this.FFPS.ffp8 ) *D*meanMotionRad;
sstDesired(6,:)=                                                      - this.FFPS.ffp5* sin( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 ) *sqrt(2*SSC/A)*meanMotionRad  + this.FFPS.ffp1 /sqrt(2*SSC*A)*   cos( sqrt(2*SSC/A)*meanMotionRad* time + this.FFPS.ffp7 ) *sqrt(2*SSC/A)*meanMotionRad;

%%^^^^^^ Re-implement formation mode selection...
% Implement based on switches
% Too confusing the way it is right now
this.StateDesired = sstDesired;
% if time<50
%   fprintf('\n %f %f %f %f %f %f',this.StateDesired);
% end
end % Function FlightControl.updateStateDesired.
