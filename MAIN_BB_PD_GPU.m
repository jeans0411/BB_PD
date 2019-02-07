% MAIN file for Bond Based Peridynamic Analysis (BB_PD) code using bond lists
%%%%%%%%%%%%% https://github.com/mhobbs18/BB_PD %%%%%%%%%%%%%%%
clear all;close all;clc;start=tic;
% global Totalnodes;
dbstop if error;
profile on;
%% Read in input data and specify material point coordinates
tic
addpath('C:\Users\ORNL-PC\Box Sync\PD Lab 0.9\Sub_Functions');
Input_Data; 

nt=1000*1000;    
fprintf('Total Loading Period: %.2fs \n', nt*dt);
Max_Force=0.5;           % Load factor
% Number of time steps (10,000 for speed testing)
countmin=1;   % Counter for determining previous time step when restarting simulations
Input_Data_time=toc;
fprintf('Input data complete in %fs \n', Input_Data_time)

%% Determine the nodes inside the horizon of each material point, build bond lists, and determine undeformed length of every bond
tic
[nodefamily,nfpointer,UndeformedLength,NumFamMembVector,NumFamMembVector2,bondlist,bondlist2,Totalbonds]=BuildHorizons(Totalnodes,coordinates,delta);
Horizons=toc;
fprintf('Horizons complete in %fs \n', Horizons)
bondlistGPU=gpuArray(bondlist);


%% Calculate volume correction factors
tic
[fac]=VolumeCorrection(Totalbonds,UndeformedLength,delta,radij);
VolumeCorrection=toc;
fprintf('Volume correction factors complete in %fs \n', VolumeCorrection)

%% Calculate bond type and bond stiffness (plus stiffness corrections)
tic
[const,BondType,BFmultiplier]=BondTypefunc(Totalbonds,bondlist,NumFamMembVector,MaterialFlag,c_concrete,c_steel,NeighbourhoodVolume,Volume);
BondTypeandStiffness=toc;
fprintf('Bond type and stiffness complete in %fs \n', BondTypeandStiffness)

%% Time integration
tic
[fail,disp,Stretch,Displacement_node,ReactionForce,countmin,LOAD_step]=TimeIntegration_GPU(Totalbonds,Totalnodes,Nod,nt,countmin,bondlistGPU,UndeformedLength,coordinates,BondType,Critical_ts_conc,Critical_ts_steel,const,Volume,fac,bodyforce,Max_Force,BFmultiplier,damping,dens,ConstraintFlag,dt,NumFamMembVector,CrossSectionFlag,Critical_ts_steel_elastic,dx);
TimeIntegrationTimer=toc;
fprintf('Time integration complete in %fs \n', TimeIntegrationTimer)
%% Simulation timing
Total=toc(start)/60/60; 
fprintf('Simulation complete in %.2f hour \n', Total)

%% Display results
tt=nt+1;
[Bond_damage]=Damage(100,tt,NumFamMembVector,fail,Totalnodes,Totalbonds,coordinates,disp,bondlistGPU,CrossSectionFlag);   % Damage - Calculate percentage of broken peridynamic bonds for every node
Displacement(disp,coordinates,MaterialFlag);                                                   % Plot deformed shape of object under analysis
DisplacementVsTime(nt,Displacement_node,LOAD_step,dx);                                                      % Plot displacement of selected node against time
%ForceVsDisplacement(Displacement_node,ForceHistory);                                          % Plot applied force against time
%StretchDisplay(Totalnodes,Stretch,Totalbonds,bondlist,coordinates,disp);                       % Plot stretch of every bond
% ReactionForceVsDisplacement(ReactionForce,Displacement_node,LOAD_step);                                  % Plot reaction force against displacement
profile viewer
% saveas(1,'f1.fig','fig');
% saveas(2,'f2.fig','fig');
% saveas(3,'f3.fig','fig');
% saveas(4,'f4.fig','fig');
% saveas(5,'f5.fig','fig');
% saveas(6,'f6.fig','fig');
% saveas(7,'f7.fig','fig');
% 
% save result7;
% system('shutdown -s');

%% Stress and Strain
% [StrainTensor]=Strainfunc(coordinates,disp,Totalbonds,bondlist,NumFamMembVector,nodefamily,nfpointer);
% [StressTensor]=Stressfunc(Totalnodes,Nod,MaterialProperties,coordinates,disp,StrainTensor,MaterialFlag);

%% Output simulation data

%% Additional time
% [Displacement_node,nt]=Additional_Time(Displacement_node,nt,countmin); % Specify additional simulation time
