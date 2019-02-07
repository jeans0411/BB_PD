function [fail,disp,Stretch,Displacement_node,ReactionForce,countmin,LOAD_step]=TimeIntegration(Totalbonds,Totalnodes,Nod,nt,countmin,bondlist,UndeformedLength,coordinates,BondType,Critical_ts_conc,Critical_ts_steel,const,Volume,fac,bodyforce,Max_Force,BFmultiplier,damping,dens,ConstraintFlag,dt,NumFamMembVector,CrossSectionFlag,Critical_ts_steel_elastic,dx)
% Time integration using a Forward Difference (FD) Backward Difference (BD) scheme (FD_BD)
% CORE Function
%% Initialise 
% global NumFamMembVector2;
fail=gpuArray.ones(Totalbonds,1);                      
Stretch=gpuArray.zeros(Totalbonds,1);
% DeformedLength=zeros(Totalbonds,1);
disp=gpuArray.zeros(Totalnodes,Nod);                         % Displacement of a material point
disp_forward=gpuArray.zeros(Totalnodes,Nod);
v=gpuArray.zeros(Totalnodes,Nod);                            % Velocity of a material point
v_forward=gpuArray.zeros(Totalnodes,Nod);
a=gpuArray.zeros(Totalnodes,Nod);                            % Acceleration of a material point
Displacement_node=gpuArray.zeros(nt,1);                      % History of displacement at selected node
% ForceHistory=zeros(nt,1);                           % History of applied force
ReactionForce=gpuArray.zeros(nt,Nod);                          
LOAD_step=gpuArray.zeros(nt,Nod);
nt0=sum([1:nt].^0.5);
failnum=gpuArray.zeros(nt,1);
dens=gpuArray(dens);
tic;
counter=0;
bondlistGPU=gpuArray(bondlist);
coordinates=gpuArray(coordinates);
for tt=countmin:nt
    countmin=tt+1; % If a checkpoint file is loaded, the simulation needs to start on the next iteration         
    %% Calculate Bond Forces and perform time integration - Forward Difference and Backward Difference scheme       
    Nforce=gpuArray.zeros(Totalnodes,Nod);  % Nodal force - initialise for every time step
    DisplacedCoordinates=coordinates+disp;  %gpu
    nt1=sum((nt-[0:tt-1]).^0.5);
    loadfactor=0.5*(1-cos(pi*nt1/nt0));
    bodyforce2=bodyforce*loadfactor;
%     [DeformedLength,Xdeformed,Ydeformed,Zdeformed,Stretch]=DeformedLengthfunc(Totalbonds,bondlist,UndeformedLength,DeformedLength,coordinates,disp,Xdeformed,Ydeformed,Zdeformed);
    [Nforce,fail,LOAD]=BondForcesGPU(Nforce,fail,BondType,Critical_ts_conc,Critical_ts_steel,const,Volume,fac,...
        bodyforce2,Max_Force,bondlistGPU,BFmultiplier,Critical_ts_steel_elastic,DisplacedCoordinates,UndeformedLength);  %all gpu
%     Nforce=gather(Nforce);
    LOAD_step(tt,:)=LOAD;   
    a(:,:)=(Nforce-damping*v(:,:))./((dens(:,:)*[1 1 1]));        % Acceleration for time:-   tt
    a(ConstraintFlag==0)=0;                                % Apply constraints
    v_forward(:,:)=v(:,:)+(a(:,:)*dt);                     % Velocity for time:-       tt + 1dt
    disp_forward(:,:)=disp(:,:) + (v_forward(:,:)*dt);     % Displacement for time:-   tt + 1dt
    v(:,:)=v_forward(:,:);                                 % Update
    disp(:,:)=disp_forward(:,:);                           % Update
    
    % Calculating the percentage of progress of time integration
    perc_progress=round((tt/nt)*100,5);


%     if perc_progress==1
%         onepercent=toc;
%         fprintf('1 percent of time integration complete in %fs \n', onepercent)
%         pause
%     end

    %% Save results
    
    
%     if mod(tt,1)==0
counter=counter+1;
Displacement_node(counter,1)=disp(14*15+8,1);          % Save displacement of defined node for plot of Displacement vs Time

% [StrainTensor]=Strainfunc(coordinates,disp,Totalbonds,bondlist,NumFamMembVector,nodefamily,nfpointer);
% [StressTensor]=Stressfunc(Totalnodes,Nod,MaterialProperties,coordinates,disp,StrainTensor,MaterialFlag);
% 
% StressHistory(counter,1)=StressTensor(66321,1,1);
%     end


% ForceHistory(tt,1)=Force;                   % Force
NforceTemp=gather(Nforce);
NforceTemp(ConstraintFlag==1)=0;
ReactionForce(tt,1)=sum(NforceTemp(:,1));
ReactionForce(tt,2)=norm(NforceTemp(:,2));
ReactionForce(tt,3)=norm(NforceTemp(:,3));

% Plot displacement of selected node against time
%ForceVsDisplacement(Displacement_node,ForceHistory);                                          % Plot applied force against time
%StretchDisplay(Totalnodes,Stretch,Totalbonds,bondlist,coordinates,disp);                       % Plot stretch of every bond
failnum(tt)=Totalbonds-sum(fail);

tt0=0;
if mod(round(perc_progress,5),1)==0
    tt0=tt0+1;
    [Bond_damage]=Damage(perc_progress,tt0,NumFamMembVector,fail,Totalnodes,Totalbonds,coordinates,disp,bondlistGPU,CrossSectionFlag);   % Damage - Calculate percentage of broken peridynamic bonds for every node
    ReactionForceVsDisplacement(ReactionForce,Displacement_node,LOAD_step,dx,failnum);
    if size(find(Bond_damage>0.99),1)>2000
        return;
    end
    %     Displacement(disp,coordinates,MaterialFlag);                                                   % Plot deformed shape of object under analysis
    %     DisplacementVsTime(nt,Displacement_node)
end
% if mod(tt,2500)==0
%     save(['D:\PhD\2 Code\BB_PD\Output\Workspace_snapshot_',num2str(tt),'.mat']); % Save workspace to local computer every 2500 time steps
% end

if mod(perc_progress,0.025)==0
    toc2=toc;tic;
    fprintf('Completed %.2f%% of time integration in %.3fs \n', perc_progress,toc2)
end

end

end
