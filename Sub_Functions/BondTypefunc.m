
function [const,BondType,BFmultiplier]=BondTypefunc(Totalbonds,bondlist,NumFamMembVector,MaterialFlag,c_concrete,c_steel,NeighbourhoodVolume,Volume)

% Determine the bond type and stiffness connecting node pairs (e.g. is it a concrete or
% steel bond?)

const=zeros(Totalbonds,1);              % bond const, Initialise vector
BondType=zeros(Totalbonds,1);       % 0,1,2; c2c? c2s, s2s
BFmultiplier=zeros(Totalbonds,1);   % 1,3; c2c s2s; c2s


for i=1:Totalbonds
        
    nodei=bondlist(i,1);
    nodej=bondlist(i,2);
        
 
    % Concrete to concrete - if materialflag i=0 and j=0 then c=concrete
    if MaterialFlag(nodei,1)==0 && MaterialFlag(nodej,1)==0
        ctemp=c_concrete;
        BondType(i)=0;
    % Concrete to steel - if materialflag i=0 and j=1 then c=concrete   
    elseif MaterialFlag(nodei,1)==0 && MaterialFlag(nodej,1)==1
        ctemp=c_concrete;
        BondType(i)=1;
    % Concrete to steel - if materialflag i=1 and j=0 then c=concrete     
    elseif MaterialFlag(nodei,1)==1 && MaterialFlag(nodej,1)==0
        ctemp=c_concrete;
        BondType(i)=1;
    % Steel to steel - if materialflag i=1 and j=1 the c=steel   
    elseif MaterialFlag(nodei,1)==1 && MaterialFlag(nodej,1)==1
        ctemp=c_steel;
        BondType(i)=2;
    end


    % Calculate stiffening factor - surface corrections for 3D problem
    NeighbourhoodVolume_i=NumFamMembVector(nodei)*Volume;                                             % Neighbourhood area for Node 'i'
    NeighbourhoodVolume_cnode=NumFamMembVector(nodej)*Volume;                                         % Neighbourhood area for Node 'j'
    StiffeningFactor=(2*NeighbourhoodVolume)/(NeighbourhoodVolume_i+NeighbourhoodVolume_cnode);
    const(i)=StiffeningFactor*ctemp;                                                                      % Correct the bond stiffness

end

% Bond force multiplier
BFmultiplier(BondType==1)=1.5;
BFmultiplier(BondType==0 | BondType==2)=1;

end