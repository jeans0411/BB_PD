% Input Data

%% Dimensions

Nod=3;          % Number of degrees of freedom (1,2,3)

Length=1;       % Length (m) - x
Width=0.1;      % Width (m) - y
Height=0.2;     % Height (m) - z

Ndiv_x=150;     % Number of divisions in x-direction
Ndiv_y=15;      % Number of divisions in y-direction
Ndiv_z=30;      % Number of divisions in z-direction

dx=Length/Ndiv_x;   % Spacing between material points in x-direction
dy=dx;              % Spacing between material points in y-direction
dz=dx;              % Spacing between material points in y-direction
% Width=Ndiv_y/Ndiv_x*Length;      % Width (m) - y
% Height=Ndiv_z/Ndiv_x*Length;     % Height (m) - z

Totalnodes=Ndiv_x*Ndiv_y*Ndiv_z;    % Total number of nodes
[coordinates]=MaterialPointCoordinates(Totalnodes,Nod,Ndiv_y,Ndiv_x,Ndiv_z,dx,dy,dz); % Define material point coordinates for main body
% loop x, then y, then z 

%% Steel rebar
MaterialFlag=zeros(Totalnodes,1);       % Create flag to identify steel and concrete nodes Concrete=0 Steel=1

% Bar 1
for i=(Ndiv_x*(14*15+7)+1):Ndiv_x*(14*15+8)
   MaterialFlag(i,1)=1; 
end
for i=(Ndiv_x*(15*15+7)+1):Ndiv_x*(15*15+8)
   MaterialFlag(i,1)=1; 
end
% for i=((Ndiv_x*378)+1):(Ndiv_x*380)
%    MaterialFlag(i,1)=1; 
% end
% for i=((Ndiv_x*393)+1):(Ndiv_x*395)
%    MaterialFlag(i,1)=1; 
% end

% Bar 2
% for i=((Ndiv_x*370)+1):(Ndiv_x*372)
%    MaterialFlag(i,1)=1; 
% end
% for i=((Ndiv_x*385)+1):(Ndiv_x*387)
%    MaterialFlag(i,1)=1; 
% end
% for i=((Ndiv_x*400)+1):(Ndiv_x*402)
%    MaterialFlag(i,1)=1; 
% end

%% Define loading plates, supports, constraints, loads (boundary conditions)
% Loading Plate
% Supports
% Constraint flag
ConstraintFlag=ones(Totalnodes,Nod);
% cantliver at fixed end
% for i=1:Ndiv_x:(Totalnodes-Ndiv_x+1)
%   ConstraintFlag(i,1)=0;
% end
% for i=2:Ndiv_x:(Totalnodes-Ndiv_x+2)
%   ConstraintFlag(i,1)=0;
% end
% for i=3:Ndiv_x:(Totalnodes-Ndiv_x+3)
%   ConstraintFlag(i,1)=0;
% end

for i=(Ndiv_x*(14*15+7)+1):Ndiv_x*(14*15+7)+1
   ConstraintFlag(i,:)=0; 
end
for i=(Ndiv_x*(15*15+7)+1):Ndiv_x*(15*15+7)+1
   ConstraintFlag(i,:)=0; 
end

for i=(Ndiv_x*(14*15+8)-0):Ndiv_x*(14*15+8)-0
   ConstraintFlag(i,2:3)=0;
end
for i=(Ndiv_x*(15*15+8)-0):Ndiv_x*(15*15+8)-0
   ConstraintFlag(i,2:3)=0; 
end


% for i=1:Ndiv_x:(Ndiv_x*Ndiv_y-1)
%   ConstraintFlag(i,3)=0;
% end
% 
% for i=1:Ndiv_x*Ndiv_x:(Totalnodes)
%   ConstraintFlag(i,2)=0;
% end

Max_Force=1.6;           % Load factor
Build_up=0;                 % Build load up over defined number of time steps
damping=4e6*5;              % Damping coefficient

% Application of force
bodyforce=zeros(Totalnodes,Nod);
% cantliver at load end
% for i=Ndiv_x:Ndiv_x:Totalnodes
%     bodyforce(i,1)=1e9/3;
% end
% for i=Ndiv_x-1:Ndiv_x:Totalnodes
%     bodyforce(i,1)=1e9/3;
% end
% for i=Ndiv_x-2:Ndiv_x:Totalnodes
%     bodyforce(i,1)=1e9/3;
% end

for i=(Ndiv_x*(14*15+8)-8):Ndiv_x*(14*15+8)-0
   bodyforce(i,1)=3.078e10/4; 
end
for i=(Ndiv_x*(15*15+8)-8):Ndiv_x*(15*15+8)-0
   bodyforce(i,1)=3.078e10/4; 
end

  

%% Material properties, and peridynamic parameters

MaterialProperties.E_concrete=22e9;                                % Young's modulus
MaterialProperties.E_steel=200e9;                                  % Young's modulus

MaterialProperties.v_concrete=0.2;                                 % Poisson's ratio
MaterialProperties.v_steel=0.3;                                    % Poisson's ratio

MaterialProperties.G_concrete=8.8e9;                               % Shear modulus
MaterialProperties.G_steel=78e9;                                   % Shear modulus
                                                   
MaterialProperties.EffectiveModulusConcrete = MaterialProperties.E_concrete/((1-2*MaterialProperties.v_concrete)*(1+MaterialProperties.v_concrete));
MaterialProperties.EffectiveModulusSteel = MaterialProperties.E_steel/((1-2*MaterialProperties.v_steel)*(1+MaterialProperties.v_steel));

MaterialProperties.Dens_concrete=2400*1e3;            % Density concrete (kg/m^3)
MaterialProperties.Dens_steel=7850*1e3;               % Density steel (kg/m^3)

dens=zeros(Totalnodes,1);
for i=1:Totalnodes
    if MaterialFlag(i,:)==0
    % Concrete
    dens(i,1)=MaterialProperties.Dens_concrete;
    elseif MaterialFlag(i,:)==1
    % Steel
    dens(i,1)=MaterialProperties.Dens_steel;        
    end
end

delta=2.1*dx;                                % Horizon, pi
NeighbourhoodVolume=(4*pi*delta^3)/3;       % Neighbourhood area/volume for node contained within material bulk
Volume=dx^3;                                % Cell volume/area
radij=dx/2;                                 % Material point radius

c_concrete=(12*MaterialProperties.E_concrete)/(pi*delta^4);        % Bond stiffness
c_steel=(12*MaterialProperties.E_steel)/(pi*delta^4);              % Bond stiffness  

Critical_ts_conc=0.000233;                  % Critical tensile stretch - concrete
Critical_ts_steel=0.03;                     % critical tensile stretch - steel   
Critical_ts_steel_elastic=60/29000;   % critical elastic tensile stretch - steel  

CrossSectionFlag=zeros(Totalnodes,1);

for i=1:Totalnodes
    if coordinates(i,2)==(1/150)*8   %select a X-Z cross section
        CrossSectionFlag(i,1)=1; % Identify nodes located in cross-section
    end        
end
%% Time Step - Look in first year report and implement defined procedure 

SF=2;  %2                                                             % Time step safety factor - to reduce time step size
dt=(0.8*sqrt(2*MaterialProperties.Dens_concrete*dx/(pi*delta^2*dx*c_concrete)))/SF;    % Minimum stable time step
nt=15000;                                                           % Number of time steps (10,000 for speed testing)
 
%% Other

countmin=1;   % Counter for determining previous time step when restarting simulations

%% Output input data into readable text file 
