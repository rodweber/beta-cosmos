function aerototalforcevector = getWindPressureVector(wind,panelSurface,noxpanels,noypanels,nozpanels,rollAngles,pitchAngles,yawAngles,rho,v,Tatmos)
  
	aerototalforcevector =zeros(3,size(rollAngles,2),size(pitchAngles,2),size(yawAngles,2));
 
  if norm(wind)==0
    return
  end
  rotspeed=30;draw=0;
  windPressure=norm(wind);
  wind        =wind/windPressure;
%   fprintf('computing aerodynamics');  
  %% must be in dimensions of force, i.e. N
  axislength=1.1*windPressure*panelSurface ;
  Ry90=roty(90);
  %Rz90=rotz(90);
  Rx90=rotx(90);
  Iz = [0 0 0.66*axislength]';
  Ix=Ry90*Iz;
  Iy=Rx90*Iz;
  %% zpanel
  pz1  = [axislength*0.9,axislength*0.9,0];pz2 = [axislength*0.9,-axislength*0.9,0];pz3 = [-axislength*0.9,-axislength*0.9,0];pz4 = [-axislength*0.9,axislength*0.9,0];
  pz12 = 0.33*pz1; pz22 = 0.33*pz2; pz32 = 0.33*pz3; pz42 = 0.33*pz4;
  pz13 = 0.66*pz1; pz23 = 0.66*pz2; pz33 = 0.66*pz3; pz43 = 0.66*pz4;
  thetaaero=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  phiaero=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  CD=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  CL=zeros(size(yawAngles,2),size(pitchAngles,2),size(rollAngles,2));
  aeroforcevectorz=[0 0 0]';
  aeroforcevectorx=[0 0 0]';
  aeroforcevectory=[0 0 0]';
  for k=1:size(yawAngles,2) %% yaw
    for j=1:size(pitchAngles,2) %% pitch
      for i=1:size(rollAngles,2) %% roll
        %% rotation matrizes
        RzY =[cosd(yawAngles(k)) -sind(yawAngles(k)) 0; sind(yawAngles(k)) cosd(yawAngles(k)) 0; 0 0 1]; %%yaw
        Ry =[cosd(pitchAngles(j))  0 sind(pitchAngles(j))  ; 0 1 0                          ; -sind(pitchAngles(j)) 0 cosd(pitchAngles(j))]; %% pitch
        RzR=[cosd(rollAngles(i)) -sind(rollAngles(i)) 0; sind(rollAngles(i)) cosd(rollAngles(i)) 0; 0 0 1]; %% roll
        if nozpanels %% zpanel
          Igz=RzY*Ry*RzR*Iz;
          [thetaaero(i,j,k),phiaero(i,j,k),Igz2]=thetaphi1(wind, Igz);
          %[aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aeroDragLiftCoef(thetaaero(i,j,k));
          [CD(i,j,k),CL(i,j,k)]= FlightControl.aeroDragLiftSentman(thetaaero(i,j,k),Tatmos,v,rho);     
          aeroforcevectorz=-wind           *    CD(i,j,k)*windPressure*panelSurface;
          ax=cross(wind,Igz2);
          if norm(ax)~=0
            liftvector = FlightControl.rodriguesRotation(wind,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end
          aeroforcevectorz=nozpanels*(-liftvector      *    CL(i,j,k)*windPressure*panelSurface + aeroforcevectorz);
        end 
        if noxpanels %% xpanel
          Igx=RzY*Ry*RzR*Ix;
          [thetaaero(i,j,k),phiaero(i,j,k),Igx2]=thetaphi1(wind, Igx);
          %[aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aeroDragLiftCoef(thetaaero(i,j,k));
          [CD(i,j,k),CL(i,j,k)]= FlightControl.aeroDragLiftSentman(thetaaero(i,j,k),Tatmos,v,rho);
          aeroforcevectorx=-wind                         *  CD(i,j,k)*windPressure*panelSurface;
          ax=cross(wind,Igx2);
          if norm(ax)~=0
            liftvector = FlightControl.rodriguesRotation(wind,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end
          aeroforcevectorx=noxpanels*(-liftvector *  CL(i,j,k)*windPressure*panelSurface + aeroforcevectorx);
        end
        if noypanels %% xpanel
          Igy=RzY*Ry*RzR*Iy;
          [thetaaero(i,j,k),phiaero(i,j,k),Igy2]=thetaphi1(wind, Igy);
          %[aerodragcoef(i,j,k),aeroliftcoef(i,j,k)]=aeroDragLiftCoef(thetaaero(i,j,k));
          [CD(i,j,k),CL(i,j,k)]= FlightControl.aeroDragLiftSentman(thetaaero(i,j,k),Tatmos,v,rho);         
          aeroforcevectory=-wind        *  CD(i,j,k)*windPressure*panelSurface;
          ax=cross(wind,Igy2) ;
          if norm(ax)~=0
            liftvector = FlightControl.rodriguesRotation(wind,ax,90/180*pi);
          else
            liftvector = [0 0 0]';
          end            
          aeroforcevectory=noypanels*(-liftvector  *  CL(i,j,k)*windPressure*panelSurface + aeroforcevectory);
        end
        aerototalforcevector(:,i,j,k)=aeroforcevectorz+aeroforcevectorx+aeroforcevectory;
        %%compute here the derivative TBC
        
        %%draw
        if draw
          vectarrow([0 0 0],wind*windPressure*panelSurface);hold on;text(wind(1)*windPressure*panelSurface,wind(2)*windPressure*panelSurface,wind(3)*windPressure*panelSurface,"wind",'HorizontalAlignment','left','FontSize',6);
          axis([-axislength axislength -axislength axislength -axislength axislength]);
          if nozpanels
            %%               
            pg = [(RzY*Ry*RzR*pz1')' ; (RzY*Ry*RzR*pz2')' ; (RzY*Ry*RzR*pz3')' ; (RzY*Ry*RzR*pz4')' ; (RzY*Ry*RzR*pz1')'];
            pg2 = [(RzY*Ry*RzR*pz12')' ; (RzY*Ry*RzR*pz22')' ; (RzY*Ry*RzR*pz32')' ; (RzY*Ry*RzR*pz42')' ; (RzY*Ry*RzR*pz12')'];
            pg3 = [(RzY*Ry*RzR*pz13')' ; (RzY*Ry*RzR*pz23')' ; (RzY*Ry*RzR*pz33')' ; (RzY*Ry*RzR*pz43')' ; (RzY*Ry*RzR*pz13')'];
            %vectarrow([0 0 0],totalforcevectorz(:,i,j,k));hold on;
            vectarrow([0 0 0],aeroforcevectorz);hold on;text(aeroforcevectorz(1),aeroforcevectorz(2),aeroforcevectorz(3),"aeroforce",'HorizontalAlignment','left','FontSize',6);
            vectarrow([0 0 0],Igz);hold on;text(Igz(1),Igz(2),Igz(3),"normal",'HorizontalAlignment','left','FontSize',6);
            line(pg(:,1), pg(:,2), pg(:,3));line(pg2(:,1), pg2(:,2), pg2(:,3));line(pg3(:,1), pg3(:,2), pg3(:,3));hold on;
            axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          if noxpanels
            %% panel and normal
            pgx = [(RzY*Ry*RzY*-Ry90*pz1')' ; (RzY*Ry*RzR*-Ry90*pz2')' ; (RzY*Ry*RzR*-Ry90*pz3')' ; (RzY*Ry*RzR*-Ry90*pz4')' ; (RzY*Ry*RzR*-Ry90*pz1')'];
            pgx2 = [(RzY*Ry*RzR*-Ry90*pz12')' ; (RzY*Ry*RzR*-Ry90*pz22')' ; (RzY*Ry*RzR*-Ry90*pz32')' ; (RzY*Ry*RzR*-Ry90*pz42')' ; (RzY*Ry*RzR*-Ry90*pz12')'];
            pgx3 = [(RzY*Ry*RzR*-Ry90*pz13')' ; (RzY*Ry*RzR*-Ry90*pz23')' ; (RzY*Ry*RzR*-Ry90*pz33')' ; (RzY*Ry*RzR*-Ry90*pz43')' ; (RzY*Ry*RzR*-Ry90*pz13')'];
            %vectarrow([0 0 0],totalforcevector(:,i,j,k)); hold on;
            vectarrow([0 0 0],Igx);hold on;text(Igx(1),Igx(2),Igx(3),"normalx",'HorizontalAlignment','left','FontSize',6);hold on;
            vectarrow([0 0 0],aeroforcevectorx);hold on;text(aeroforcevectorx(1),aeroforcevectorx(2),aeroforcevectorx(3),"aeroforcex",'HorizontalAlignment','left','FontSize',6);
            line(pgx(:,1), pgx(:,2), pgx(:,3));line(pgx2(:,1), pgx2(:,2), pgx2(:,3));line(pgx3(:,1), pgx3(:,2), pgx3(:,3));hold on;
            axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          if noypanels
            %% panel and normal
            pgy = [(RzY*Ry*RzR*-Rx90*pz1')' ; (RzY*Ry*RzR*-Rx90*pz2')' ; (RzY*Ry*RzR*-Rx90*pz3')' ; (RzY*Ry*RzR*-Rx90*pz4')' ; (RzY*Ry*RzY*-Rx90*pz1')'];
            pgy2 = [(RzY*Ry*RzR*-Rx90*pz12')' ; (RzY*Ry*RzR*-Rx90*pz22')' ; (RzY*Ry*RzR*-Rx90*pz32')' ; (RzY*Ry*RzR*-Rx90*pz42')' ; (RzY*Ry*RzR*-Rx90*pz12')'];
            pgy3 = [(RzY*Ry*RzR*-Rx90*pz13')' ; (RzY*Ry*RzR*-Rx90*pz23')' ; (RzY*Ry*RzR*-Rx90*pz33')' ; (RzY*Ry*RzR*-Rx90*pz43')' ; (RzY*Ry*RzR*-Rx90*pz13')'];
            %vectarrow([0 0 0],totalforcevector(:,i,j,k));hold on;                  
            vectarrow([0 0 0],Igy);hold on;text(Igy(1),Igy(2),Igy(3),"normaly",'HorizontalAlignment','left','FontSize',6);
            vectarrow([0 0 0],aeroforcevectory);hold on;text(aeroforcevectory(1),aeroforcevectory(2),aeroforcevectory(3),"aeroforcey",'HorizontalAlignment','left','FontSize',6);
            line(pgy(:,1), pgy(:,2), pgy(:,3));line(pgy2(:,1), pgy2(:,2), pgy2(:,3));line(pgy3(:,1), pgy3(:,2), pgy3(:,3));hold on;                        axis([-axislength axislength -axislength axislength -axislength axislength]);
          end
          %text(0,0,0,strcat('scalingfactor: ',int2str(scalingfactor)),'HorizontalAlignment','left','FontSize',6);
          xlabel('fx [N]');ylabel('fy [N]');zlabel('fz [N]');
          axis([-axislength axislength -axislength axislength -axislength axislength]);
          hold off;
          pause(1/rotspeed);
        end
      end
      %input('qaz')
    end
  end
%   fprintf(' - done\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                       
function [theta,phi,Ig2]=thetaphi1(refvec, vec)
%% this function computes theta,phi, Ig2
  theta = atan2d(norm(cross(refvec,vec)), dot(refvec,vec));                                                                                                       
  if theta>90                                                                                                                                                      
    theta=180-theta;                                                                                                                                             
    Ig2=-vec;                                                                                                                                                     
  else                                                                                                                                                             
    Ig2=vec;                                                                                                                                                     
  end                                                                                                                                                              
  phi=atand( (refvec(3)-vec(3)) / (refvec(2)-vec(2)+1e-30) );  
end     

