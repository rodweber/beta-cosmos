function [forceVector,rollAngleOpt,pitchAngleOpt,yawAngleOpt]=findBestAttitude(totalForceVector,u,rollAngles,pitchAngles,...
                                                 yawAngles,oldRollAngleOpt,oldPitchAngleOpt,oldYawAngleOpt) 
%findBestAerodynamicAngles to find best angle
%%%%%%% input variables
%% totalforcevector
%% controlvector
%% rollAngles
%% pitchAngles
%% yawAngles
%% oldRollAngleOpt
%% oldPitchAngleOpt
%% oldYawAngleOpt
%%%%%%% output variables
%% forcevector
%% rollAngleOpt
%% pitchAngleOpt
%% yawAngleOpt

%u=[-3e-6 -15e-6 -2e-5]'

  thetaRange=90;
  %% is this a universally good thetaRange?
  %% What is the relation to angle granularity?
  %% Do we need a formula like 1.5xgranularity,
  %% or a check that thetaRange is larger than granularity?
  
  theta=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  phi=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));

  l=1;
  
  for k=1:size(yawAngles,2) %% yaw
    for j=1:size(pitchAngles,2) %% pitch
      for i=1:size(rollAngles,2) %% roll
        %rollAngles(i)
        %pitchAngles(j)
        %yawAngles(k)
        %norm(totalForceVector(:,i,j,k))
        %input('XXX')
        if norm(totalForceVector(:,i,j,k))~=0 
          %[theta(i,j,k),phi(i,j,k)]=thetaphi(totalForceVector(:,i,j,k)/norm(totalForceVector(:,i,j,k)),u/norm(u));        
          theta(i,j,k) = atan2d(norm(cross(totalForceVector(:,i,j,k),u)), dot(totalForceVector(:,i,j,k),u));
          %theta(i,j,k)
          %u
          if theta(i,j,k)<=thetaRange %% find all force vectors that have theta<thetarange deg
              goodTheta(l)=theta(i,j,k);
              goodThetai(l)=i; goodThetaj(l)=j; goodThetak(l)=k;
              goodRollAngle(l)=rollAngles(i); goodPitchAngle(l)=pitchAngles(j); goodYawAngle(l)=yawAngles(k);
              %! some proper weighting factors are needed
              %meritFactor(l)=1000/(  (1e-5+abs(theta(i,j,k))) * (1+wrapTo360(oldRollAngleOpt-rollAngles(i))) * (1+wrapTo360(oldPitchAngleOpt-pitchAngles(j))) * (1+wrapTo360(oldYawAngleOpt-yawAngles(k))) );
              mTheta=1e-7;mRollAngle=1000;mPitchAngle=1000;mYawAngle=1000;
              %meritFactor(l)=1/(mTheta+4*abs(theta(i,j,k))) + 1/(mRollAngle+wrapTo360(oldRollAngleOpt-rollAngles(i))) +1/(mPitchAngle+2*wrapTo360(oldPitchAngleOpt-pitchAngles(j))) + 1/(mYawAngle+wrapTo360(oldYawAngleOpt-yawAngles(k)));
              meritFactorInverse(l)=theta(i,j,k);
              l=l+1;
              %vectarrow([0 0 0],u/norm(u));hold on;
              %vectarrow([0 0 0],totalForceVector(:,i,j,k)/norm(totalForceVector(:,i,j,k)));
              %hold on;axis equal;

          end
        end
        %input('xxx')
      end
    end
  end
  if not(exist('meritFactorInverse')) %% enter if meritFactor does not exist

%    rollAngleOpt=oldRollAngleOpt;
%    pitchAngleOpt=oldPitchAngleOpt;
%    yawAngleOpt=oldYawAngleOpt;
%    ii=find(rollAngles==oldRollAngleOpt);
%    jj=find(pitchAngles==oldPitchAngleOpt);
%    kk=find(yawAngles==oldYawAngleOpt);
%    forceVector=squeeze(totalForceVector(:,ii,jj,kk));
 
    rollAngleOpt=0; pitchAngleOpt=0; yawAngleOpt=0;
    forceVector=[0 0 0]';
    %u'
    %fprintf('\n findBestAerodynamicAngles: problem with theta\n')
    %input('XXX')
  elseif 0 %oldRollAngleOpt==0 && oldPitchAngleOpt==0 && oldYawAngleOpt==0
    [~,optIndex]=max(goodTheta);
    rollAngleOpt=goodRollAngle(optIndex);
    pitchAngleOpt=goodPitchAngle(optIndex);
    yawAngleOpt=goodYawAngle(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
    %fprintf('\nfindbestAerodynamicAngles - elseif\n');
  else
    [~,optIndex]=min(meritFactorInverse);
    %meritFactor
    %a=min(meritFactorInverse)
    %idx=find(meritFactorInverse==a)
    %goodTheta
    rollAngleOpt=goodRollAngle(optIndex);
    pitchAngleOpt=goodPitchAngle(optIndex);
    yawAngleOpt=goodYawAngle(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
    %fprintf('\nfindbestAerodynamicAngles - else\n');
      %figure
      %histogram(meritFactor)

  end  
  %!compute force for uncertain angles
  uncertainty=0;
  uncertaintyX=uncertainty;   uncertaintyY=uncertainty;  uncertaintyZ=uncertainty;
  variation=[(1+uncertaintyX*(1-2*rand)); (1+uncertaintyY*(1-2*rand)); (1+uncertaintyZ*(1-2*rand))];
  forceVector=forceVector.*variation;
  %hold off;
    %magcv=sqrt(controlvector(1)^2+controlvector(2)^2+controlvector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magcv,controlvector(1)/magcv,controlvector(2)/magcv,controlvector(3)/magcv);
    %magfv=sqrt(forcevector(1)^2+forcevector(2)^2+forcevector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magfv,forcevector(1)/magfv,forcevector(2)/magfv,forcevector(3)/magfv);
    %fprintf('-------')

%{
    figure
      histogram(meritfactor)
    figure
      plot(meritfactor)
    figure
      subplot(2,1,1)
        histogram(theta)
      subplot(2,1,2)
        histogram(phi)
    figure
      for n=1:l-1
        vectarrow([0 0 0],totalforcevector(:,squeeze(goodthetai(n)),squeeze(goodthetaj(n)),squeeze(goodThetak(n))));
        hold on;axis equal;
      end
      hold off;
    Iz=[0 0 1]';
    figure
      for n=1:l-1
        RzY =[cosd(goodyawAngle(n)) -sind(goodyawAngle(n)) 0; sind(goodyawAngle(n)) cosd(goodyawAngle(n)) 0; 0 0 1]; %%roll
        Ry =[cosd(goodpitchAngle(n))  0 sind(goodpitchAngle(n))  ; 0 1 0                          ; -sind(goodpitchAngle(n)) 0 cosd(goodpitchAngle(n))]; %% pitch
        RzR=[cosd(goodrollAngle(n)) -sind(goodrollAngle(n)) 0; sind(goodrollAngle(n)) cosd(goodrollAngle(n)) 0; 0 0 1]; %% yaw
        Ig=RzY*Ry*RzR*Iz;
        vectarrow([0 0 0],Ig);
        hold on;axis equal;
      end
    hold off;
    
%}
%forceVector
%input('test1'); 
end



%{
 OLD FUNCTION


function [forceVector,rollAngleOpt,pitchAngleOpt,yawAngleOpt]=findBestAttitude(totalForceVector,controlVector,rollAngles,pitchAngles,...
                                                 yawAngles,oldRollAngleOpt,oldPitchAngleOpt,oldYawAngleOpt)
%findBestAerodynamicAngles to find best angle
%%%%%%% input variables
%% totalforcevector
%% controlvector
%% rollAngles
%% pitchAngles
%% yawAngles
%% oldRollAngleOpt
%% oldPitchAngleOpt
%% oldYawAngleOpt
%%%%%%% output variables
%% forcevector
%% rollAngleOpt
%% pitchAngleOpt
%% yawAngleOpt

  thetaRange=45;
  %% is this a universally good thetaRange?
  %% What is the relation to angle granularity?
  %% Do we need a formula like 1.5xgranularity,
  %% or a check that thetaRange is larger than granularity?
  
  theta=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  phi=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  l=1;
  for k=1:size(yawAngles,2) %% yaw
    for j=1:size(pitchAngles,2) %% pitch
      for i=1:size(rollAngles,2) %% roll
        
        %totalForceVector(:,i,j,k)'/norm(totalForceVector(:,i,j,k))
        %controlVector'/norm(controlVector)
        [theta(i,j,k),phi(i,j,k)]=thetaphi(totalForceVector(:,i,j,k)/norm(totalForceVector(:,i,j,k)),controlVector/norm(controlVector));        
        %theta(i,j,k)
        %phi(i,j,k)
        if theta(i,j,k)<thetaRange %% find all force vectors that have theta<thetarange deg
          %fprintf('-x-\n')
          goodTheta(l)=theta(i,j,k);
          goodThetai(l)=i;
          goodThetaj(l)=j;
          goodThetak(l)=k;
          goodRollAngle(l)=rollAngles(i);
          goodPitchAngle(l)=pitchAngles(j);
          goodYawAngle(l)=yawAngles(k);
          %! some proper weighting factors are needed
          meritFactor(l)=1000/(  (1e-5+abs(theta(i,j,k)))  *    (1+wrapTo360(oldRollAngleOpt-rollAngles(i)))  *(1+wrapTo360(oldPitchAngleOpt-pitchAngles(j))) * (1+wrapTo360(oldYawAngleOpt-yawAngles(k))) );
          l=l+1;
        end
        %input('xxx')
      end
    end
  end
  if not(exist('meritFactor'))
    rollAngleOpt=oldRollAngleOpt;
    pitchAngleOpt=oldPitchAngleOpt;
    yawAngleOpt=oldYawAngleOpt;
    ii=find(rollAngles==oldRollAngleOpt);
    jj=find(pitchAngles==oldPitchAngleOpt);
    kk=find(yawAngles==oldYawAngleOpt);
    ii=1;jj=1;kk=1;
    forceVector=squeeze(totalForceVector(:,ii,jj,kk));
    %controlVector/norm(controlVector)
    %forceVector/norm(forceVector)
    %fprintf('\n findBestAerodynamicAngles: problem with theta\n')
    %input('')
  elseif oldRollAngleOpt==0 && oldPitchAngleOpt==0 && oldYawAngleOpt==0
    [~,optIndex]=min(goodTheta);
    rollAngleOpt=goodRollAngle(optIndex);
    pitchAngleOpt=goodPitchAngle(optIndex);
    yawAngleOpt=goodYawAngle(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
  else
    [~,optIndex]=max(meritFactor);
    rollAngleOpt=goodRollAngle(optIndex);
    pitchAngleOpt=goodPitchAngle(optIndex);
    yawAngleOpt=goodYawAngle(optIndex);
    forceVector=squeeze(totalForceVector(:,goodThetai(optIndex),goodThetaj(optIndex),goodThetak(optIndex)));
  end  

    %magcv=sqrt(controlvector(1)^2+controlvector(2)^2+controlvector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magcv,controlvector(1)/magcv,controlvector(2)/magcv,controlvector(3)/magcv);
    %magfv=sqrt(forcevector(1)^2+forcevector(2)^2+forcevector(3)^2);
    %fprintf('\n %1.1e %1.3f %1.3f %1.3f',magfv,forcevector(1)/magfv,forcevector(2)/magfv,forcevector(3)/magfv);
    %fprintf('-------')

%{
    figure
      histogram(meritfactor)
    figure
      plot(meritfactor)
    figure
      subplot(2,1,1)
        histogram(theta)
      subplot(2,1,2)
        histogram(phi)
    figure
      for n=1:l-1
        vectarrow([0 0 0],totalforcevector(:,squeeze(goodthetai(n)),squeeze(goodthetaj(n)),squeeze(goodThetak(n))));
        hold on;axis equal;
      end
      hold off;
    Iz=[0 0 1]';
    figure
      for n=1:l-1
        RzY =[cosd(goodyawAngle(n)) -sind(goodyawAngle(n)) 0; sind(goodyawAngle(n)) cosd(goodyawAngle(n)) 0; 0 0 1]; %%roll
        Ry =[cosd(goodpitchAngle(n))  0 sind(goodpitchAngle(n))  ; 0 1 0                          ; -sind(goodpitchAngle(n)) 0 cosd(goodpitchAngle(n))]; %% pitch
        RzR=[cosd(goodrollAngle(n)) -sind(goodrollAngle(n)) 0; sind(goodrollAngle(n)) cosd(goodrollAngle(n)) 0; 0 0 1]; %% yaw
        Ig=RzY*Ry*RzR*Iz;
        vectarrow([0 0 0],Ig);
        hold on;axis equal;
      end
    hold off;
    
%}
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,phi]=thetaphi(refvec, vec)
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)) );
end



%}