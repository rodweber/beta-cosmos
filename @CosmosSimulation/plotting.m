function plotting(this, ns, meanMotionRad)
%function plotting(angles,sst,refPosChange,time,ns,meanMotion,u,e)

plotExperimentTime=1;
plot3D=1;
generalVariablesPlot=1;
controlVariablesPlot=1;

% Get path to telemetry files from CosmosSimulation object.
tmFolderPath = this.TelemetryPath;

antennaConeHalfAngle=60; %%[deg]
%plannedExperimentTimes=readmatrix('2times20minperdayfor7days.csv');
%plannedExperimentTimes=readmatrix('2times40minperdayfor7days.csv');
%plannedExperimentTimes=readmatrix('2times240minperdayfor7days.csv');
%plannedExperimentTimes=readmatrix('everyday6times30minperdayfor7daysTBC.csv');
plannedExperimentTimes=readmatrix('everyday6times30minperdayfor7daysTBC.csv');


%% read data from telemetry files
for i=1:ns
  fileControlVectorTM = strcat(tmFolderPath,filesep,'ControlVectorTM',num2str(i),'.csv');
  fileForceVectorTM = strcat(tmFolderPath,filesep,'ForceVectorTM',num2str(i),'.csv');
  fileSatPositionTM = strcat(tmFolderPath,filesep,'SatPositionTM',num2str(i),'.csv');
  fileTimeVectorTM = strcat(tmFolderPath,filesep,'TimeVectorTM',num2str(i),'.csv');
  fileSatStatesTM = strcat(tmFolderPath,filesep,'SatStatesTM',num2str(i),'.csv');
  
  tempTime=readmatrix(fileTimeVectorTM);
  tempRefPosChange=readmatrix(fileSatPositionTM);
  tempSatStates=readmatrix(fileSatStatesTM);
  tempControlVector=readmatrix(fileControlVectorTM);
  tempForceVector  =readmatrix(fileForceVectorTM);
  if i==1
    timeSteps=size(tempTime,1);
    cosmosTime=zeros(timeSteps,ns);
    sst=zeros(ns,6,timeSteps);
    angles=zeros(ns,3,timeSteps);
    refPosChange=zeros(1,3,timeSteps);
    refPosChange(1,:,:)=tempRefPosChange';
    controlVector=zeros(ns,3,timeSteps);
    forceVector=zeros(ns,3,timeSteps);
  end
  cosmosTime(:,i)=tempTime(:)';
  sst(i,:,:)=tempSatStates(:,1:6)';
  angles(i,:,:)=tempSatStates(:,7:9)';
  controlVector(i,:,:)=tempControlVector';
  forceVector(i,:,:)=tempForceVector';
end

if plotExperimentTime%% condition for experiment time valid?
  TAngle=zeros(size(cosmosTime,1),1);
  SAngle=zeros(size(cosmosTime,1),1);
  experimentTime=zeros(size(cosmosTime,1),1);
  plannedExperimentTime=zeros(size(cosmosTime,1),1);
  directionAngle=zeros(size(cosmosTime,1),1);
  

%% create function for when an experiment time is defined  
for i=1:size(cosmosTime,1)
  for j=1:size(plannedExperimentTimes,1)
    if cosmosTime(i,1)>plannedExperimentTimes(j,1) && cosmosTime(i,1)<plannedExperimentTimes(j,2)
      plannedExperimentTime(i)=1;
    end
  end
end
  
  for i=1:size(cosmosTime,1)
    %% compute angle
    directionAngle(i)= atan2d( (sst(2,3,i)-sst(1,3,i)) , sqrt((sst(2,1,i)-sst(1,1,i))^2+(sst(2,2,i)-sst(1,2,i))^2) );
    %% compute angle with the local x axis (Hill coordinate system)
    %if angles(1,2,i)<90 && angles(1,3,i)<90 && angles(2,2,i)<90 && angles(2,3,i)<90
    TPitchYawAngle=atan2d(  1,sqrt(1/tand(angles(1,2,i))^2+1/tand(angles(1,3,i))^2)  );
    SPitchYawAngle=atan2d(  1,sqrt(1/tand(angles(2,2,i))^2+1/tand(angles(2,3,i))^2)  );
    %else
    %  TpolarAngle=90;SpolarAngle=90;
    %end
    %% for each satellite, compute angle towards the other satellite
    TAngle(i)        =TPitchYawAngle-directionAngle(i);
    SAngle(i)        =SPitchYawAngle+directionAngle(i);
    %% if for each satellite, the other satellite is in view, then an experiment can be carried out
    if abs(TAngle(i))<antennaConeHalfAngle && abs(SAngle(i))<antennaConeHalfAngle
      experimentTime(i)=1;
    end
  end
  figure
  plot(cosmosTime(:,1)/2/pi*meanMotionRad,TAngle,  cosmosTime(:,1)/2/pi*meanMotionRad,SAngle,  cosmosTime(:,1)/2/pi*meanMotionRad,directionAngle);hold on;
  axis([-inf inf -10 190])
  yticks([0 45 90 135 180])
  yyaxis right;
  plot(cosmosTime(:,1)/2/pi*meanMotionRad,experimentTime,'-k','LineWidth',2);
  plot(cosmosTime(:,1)/2/pi*meanMotionRad,plannedExperimentTime,'-g','LineWidth',2);
  axis([0 inf -.1 1.1])
  %title('BETA: experiment time');
  legend('TAngle','SAngle','directionAngle','experiment time achieved (right)','experiment time defined (right)');
  hold off;
end





if generalVariablesPlot %% plot general variables
  figure
  
  subplot(2,3,1)%% roll
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(angles(i,1,:)));hold on
    names(i)=[{strcat('sat',int2str(i))}];
  end
  ylabel('roll angle [deg]');xlabel('no. of orbits');grid on;hold off;
  legend(names);
  axis([-inf inf -10 370])
  yticks([0 45 90 135 180 225 270 315 360])
  %title('BETA');
  
  subplot(2,3,2)%% pitch
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(angles(i,2,:)));hold on
  end
  ylabel('pitch angle [deg]');xlabel('no. of orbits');grid on;hold off;
  axis([-inf inf -10 190])
  yticks([0 45 90 135 180])
  
  subplot(2,3,3)%%yaw
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(angles(i,3,:)));hold on
  end
  ylabel('yaw angle [deg]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -10 370])
  yticks([0 45 90 135 180 225 270 315 360])
  
  subplot(2,3,4)%%x
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(sst(i,1,:)));hold on
  end
  ylabel('x [m]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  

  subplot(2,3,5)%%y
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(sst(i,2,:)));hold on
  end
  ylabel('y [m]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  
  
  subplot(2,3,6)%%z
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(sst(i,3,:)));hold on
  end
  ylabel('z [m]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  

end






 %% plot control and force
if controlVariablesPlot
  figure
  subplot(2,3,1)%% control vector x
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(controlVector(i,1,:)));hold on
  end
  ylabel('control vector x [?]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  
  %title('BETA');
  
  subplot(2,3,2)%% control vector y
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(controlVector(i,2,:)));hold on
  end
  ylabel('control vector y [?]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  

  subplot(2,3,3)%% control vector z
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(controlVector(i,3,:)));hold on
  end
  ylabel('control vector z [?]');xlabel('no. of orbits');grid on;hold off
  axis([-inf inf -inf inf])  

  subplot(2,3,4)%% force vector x direction
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(forceVector(i,1,:)));hold on
    names(i)=[{strcat('sat',int2str(i))}];
  end
  ylabel('force vector x [?]');xlabel('no. of orbits');grid on;hold off;
  axis([-inf inf -inf inf])  
  legend(names);
  
  subplot(2,3,5)%% force vector y direction
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(forceVector(i,2,:)));hold on
  end
  ylabel('force vector y [?]');xlabel('no. of orbits');grid on;hold off;
  axis([-inf inf -inf inf])  
  
  subplot(2,3,6)%% force vector z
  for i=1:ns
    plot(squeeze(cosmosTime(:,i)/2/pi*meanMotionRad),squeeze(forceVector(i,3,:)));hold on
  end
  ylabel('force vector z [?]');xlabel('no. of orbits');grid on;hold off;
  axis([-inf inf -inf inf])  
end

%% 3d plotting
if plot3D
  figure
  fontSize=40;
  lineWidth=1;
  for i=1:ns
    plot3(squeeze(sst(i,1,:)),squeeze(sst(i,2,:)),squeeze(sst(i,3,:)),'-','LineWidth',lineWidth);hold on;
    names(i)=[{strcat('sat',int2str(i))}];
    %plot3(xzyplane(1,:),xzyplane(2,:),xzyplane(3,:),squeeze(sst(1,i,:)),squeeze(sst(2,i,:)),squeeze(sst(3,i,:)),'LineWidth',lineWidth);hold on
  end
  xlabel('X [m]','FontSize', fontSize);ylabel('Y [m]','FontSize', fontSize);zlabel('Z [m]','FontSize', fontSize);
  %legend('this research','FontSize', fontSize);
  legend(names,'FontSize', fontSize,'units','centimeters','Position',[19 15.5 1 1])
  %legend('w/o solar pressure','w/ solar pressure','FontSize', fontSize,'units','centimeters','Position',[15 17 1 1]);
  grid on;axis equal;
  set(gcf,'units','centimeters','position',[0,0,27,20]);
  set(gca,'units','centimeters','position',[5.5 4 16 15]);
  set(gca,'FontSize',fontSize);
  %axis([-2300 150 -320 320 -320 320]);
  %xticks([-2000 -1000 0 ]);yticks([-200 200]);zticks([-200 200]);
  %view([1 -0.5 0.5]);
  %set(gcf, 'InvertHardCopy', 'off');
  %title('3d');
  %title('BETA');
  hold off;
  
  %savefig('sst3d.fig');
  %print('-painters','-dmeta','sst3d.emf')
  
end

%%
if 0
  figure
  for i=1:3
    subplot(3,1,i)
    for j=1:ns
      plot(cosmosTime(:,j)/2/pi*meanMotionRad,squeeze(e(j,i,:)));hold on;
    end
    legend;
    %title(strcat('BETA: error coordinate',num2str(i)))
  end
  hold off;
end

%% reference position change
if 0 
  refPosChange = squeeze(refPosChange);
  figure
  plot(cosmosTime(:,1),refPosChange(1,:),cosmosTime(:,1),refPosChange(2,:),cosmosTime(:,1),refPosChange(3,:));
  legend('x','y','z');
  %title('BETA: reference position change')
  hold off;
end

if 0 %% 2d printing
  figure
  plot(traub(:,1),traub(:,2),xzplane(1,:),xzplane(3,:),squeeze(sst(1,i,:)),squeeze(sst(3,i,:)))
  xlabel('X [m]','FontSize', fontSize);ylabel('Y [m]','FontSize', fontSize);
  grid on;axis equal;hold off;
  axis([-1400 100 -400 600]);
  box off;
  legend('Traub et al. [4] ','this research in-plane','this research in/out-of-plane','FontSize', fontSize,'Units','centimeters','Position',[16 16 1 1]);
  
  %set(gca,'FontSize',fontSize);
  %set(gcf,'units','centimeters','position',[0,0,25,20])
  %legend('w/o solar pressure','w/ solar pressure','FontSize', fontSize,'Units','centimeters','Position',[12 17 1 1]);
  set(gcf,'units','centimeters','position',[0,0,27,20]);
  set(gca,'units','centimeters','position',[5.5 3 20 15]);
  set(gca,'FontSize',fontSize);
  print('-painters','-dmeta','sst2d.emf')
  %title('BETA');
end
end



