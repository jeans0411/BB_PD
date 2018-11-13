% MAIN file for Bond Based Peridynamic Analysis (BB_PD) code using bond lists

clear all
close all
clc

start=tic;

%% Read in input data and specify material point coordinates
tic
addpath('D:\PhD\2 Code\BB_PD\Sub_Functions');
Input_Data                                     
Input_Data_time=toc;
fprintf('Input data complete in %fs \n', Input_Data_time)

%% Determine the nodes inside the horizon of each material point
tic
[nodefamily,nfpointer,UndeformedLength,NumFamMembVector,MaxNumFamMemb,bondlist,Totalbonds]=BuildHorizons(Totalnodes,coordinates,delta);
Horizons=toc;
fprintf('Horizons complete in %fs \n', Horizons)

%% Calculate volume correction factors
tic
[fac]=VolumeCorrection(Totalbonds,UndeformedLength,delta,radij);
VolumeCorrection=toc;
fprintf('Volume correction factors complete in %fs \n', VolumeCorrection)

%% Calculate bond type and bond stiffness (plus stiffness corrections)
tic
[c,BondType]=BondType(Totalbonds,bondlist,NumFamMembVector,MaterialFlag,c_concrete,c_steel,NeighbourhoodVolume,Volume);
BondTypeandStiffness=toc;
fprintf('Bond type and stiffness complete in %fs \n', BondTypeandStiffness)
%% Time integration
tic
[fail,disp,Stretch,Displacement_node,ReactionForce,countmin]=TimeIntegration(Totalnodes,Totalbonds,bondlist,Nod,countmin,nt,Build_up,Max_Force,NumFamMembVector,nodefamily,nfpointer,UndeformedLength,coordinates,BondType,Critical_ts_conc,Critical_ts_steel,c,fac,Volume,damping,dens,dt,ConstraintFlag,bodyforce);
TimeIntegration=toc;
fprintf('Time integration complete in %fs \n', TimeIntegration)

%% Simulation timing
Total=toc(start); 
fprintf('Simulation complete in %fs \n', Total)

%% Display results
[Bond_damage]=Damage(NumFamMembVector,fail,Totalnodes,coordinates,disp);   % Damage - Calculate percentage of broken peridynamic bonds for every node
Displacement(disp,coordinates,MaterialFlag);                               % Plot deformed shape of object under analysis
DisplacementVsTime(nt,Displacement_node);                                  % Plot displacement of selected node against time
%ForceVsDisplacement(Displacement_node,ForceHistory);                      % Plot applied force against time
Stretch(Stretch,coordinates,disp);                                              % Plot stretch of every bond
ReactionForceVsDisplacement(ReactionForce,Displacement_node);                   % Plot reaction force against displacement

%% Stress and Strain
%[StrainTensor]=Strains(coordinates,Displacement,NumFamMemb,nodefamily,Totalnodes,familypointer,maxfam,delta);
%[SressTensor]=Stress(Totalnodes,Nod,EffectiveModulusConcrete,EffectiveModulusSteel,v_concrete,v_steel,G_concrete,G_steel,coordinates,Displacement,StrainTensor,MaterialFlag);

%% Output simulation data

%% Additional time
% [Displacement_node,nt]=Additional_Time(Displacement_node,nt,countmin); % Specify additional simulation time
