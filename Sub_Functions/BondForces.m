function [Nforce,fail,LOAD]=BondForces(tt,nt0,nt,Nforce,Totalbonds,fail,BondType,Stretch,Critical_ts_conc,Critical_ts_steel,const,Volume,fac,DeformedLength,Xdeformed,Ydeformed,Zdeformed,bodyforce,Max_Force,bondlist,BFmultiplier,Critical_ts_steel_elastic)
% Calculate bond forces, make use of logical indexing

% Initialise bond force to zero for every time step
BforceX=zeros(Totalbonds,1); 
BforceY=zeros(Totalbonds,1);
BforceZ=zeros(Totalbonds,1);

fail(fail==1 & BondType==0 & Stretch>Critical_ts_conc)=0;     % Deactivate bond if stretch exceeds critical stretch   Failed = 0 
fail(fail==1 & BondType==1 & Stretch>3*Critical_ts_conc)=0;   % EMU user manual recommends that the critical stretch and bond force are multiplied by a factor of 3 for concrete to steel bonds 
fail(fail==1 & BondType==2 & Stretch>Critical_ts_steel)=0;    
% Bond remains active = 1

% Calculate X,Y,Z component of bond force
BforceX=BFmultiplier.*fail.*const.*Stretch*Volume.*fac.*(Xdeformed./DeformedLength);
BforceY=BFmultiplier.*fail.*const.*Stretch*Volume.*fac.*(Ydeformed./DeformedLength);
BforceZ=BFmultiplier.*fail.*const.*Stretch*Volume.*fac.*(Zdeformed./DeformedLength);

%check steel perfect plastic
plastic=fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic;
ultimate=fail==1 & BondType==2 & Stretch>0.01;


if isempty(single(plastic(plastic==true)))==0
    ans=find(single(plastic)==1);
    BforceX(ans)=BforceX(ans)*Critical_ts_steel_elastic./Stretch(ans);
    BforceY(ans)=BforceY(ans)*Critical_ts_steel_elastic./Stretch(ans);
    BforceZ(ans)=BforceZ(ans)*Critical_ts_steel_elastic./Stretch(ans);
end

if isempty(single(ultimate(ultimate==true)))==0
    ans=find(single(ultimate)==1);
    BforceX(ans)=BforceX(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
    BforceY(ans)=BforceY(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
    BforceZ(ans)=BforceZ(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
end



%

% BforceX(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceX(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% BforceY(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceY(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% BforceZ(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceZ(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);

% DeformedLength will be 0 intially and divison leads to 'Not-a-Number'. If Bforce = 'NaN', set to 0 
BforceX(isnan(BforceX))=0; 
BforceY(isnan(BforceY))=0;
BforceZ(isnan(BforceZ))=0;

% Calculate the nodal force for every node, iterate over the bond list
% for i=1:Totalbonds    
%     nodei=bondlist(i,1); % Node i
%     nodej=bondlist(i,2); % Node j   
%     % X-component
%     Nforce(nodei,1)=Nforce(nodei,1)+BforceX(i); % Bond force is positive on Node i 
%     Nforce(nodej,1)=Nforce(nodej,1)-BforceX(i); % Bond force is negative on Node j   
%     % Y-component
%     Nforce(nodei,2)=Nforce(nodei,2)+BforceY(i);
%     Nforce(nodej,2)=Nforce(nodej,2)-BforceY(i);
%     % Z-component
%     Nforce(nodei,3)=Nforce(nodei,3)+BforceZ(i);
%     Nforce(nodej,3)=Nforce(nodej,3)-BforceZ(i);
%         
% end

% i=1;
% while i<=Totalbonds
%     nodei=bondlist(i,1); % Node i
%     nodej=bondlist(i,2); % Node j  
%     % X-component
%     Nforce(nodei,:)=Nforce(nodei,:)+[BforceX(i) BforceY(i) BforceZ(i)]; % Bond force is positive on Node i 
%     Nforce(nodej,:)=Nforce(nodej,:)-[BforceX(i) BforceY(i) BforceZ(i)]; % Bond force is negative on Node j
% %     % Y-component
% %     Nforce(nodei,2)=Nforce(nodei,2)+BforceY(i);
% %     Nforce(nodej,2)=Nforce(nodej,2)-BforceY(i); 
% %     % Z-component
% %     Nforce(nodei,3)=Nforce(nodei,3)+BforceZ(i);
% %     Nforce(nodej,3)=Nforce(nodej,3)-BforceZ(i);  
%     i=i+1;
% end

i=1;
while i<=Totalbonds
    nodei=bondlist(i,1); % Node i
    nodej=bondlist(i,2); % Node j  
    % X-component
    Nforce(nodei,1)=Nforce(nodei,1)+BforceX(i); % Bond force is positive on Node i 
    Nforce(nodej,1)=Nforce(nodej,1)-BforceX(i); % Bond force is negative on Node j
    % Y-component
    Nforce(nodei,2)=Nforce(nodei,2)+BforceY(i);
    Nforce(nodej,2)=Nforce(nodej,2)-BforceY(i); 
    % Z-component
    Nforce(nodei,3)=Nforce(nodei,3)+BforceZ(i);
    Nforce(nodej,3)=Nforce(nodej,3)-BforceZ(i);  
    i=i+1;
end
                                          
% Add body force

nt1=sum((nt-[1:tt]+1).^0.5);

loadfactor=0.5*(1-cos(pi*nt1/nt0));

bodyforce2=bodyforce*loadfactor;
Nforce(:,:)=Nforce(:,:)+(bodyforce2(:,:)*Max_Force);
LOAD=sum((bodyforce2(:,:)*Max_Force));
end
