function [Nforce,fail,LOAD]=BondForcesGPU(Nforce,fail,BondType,Critical_ts_conc,Critical_ts_steel,const,Volume,fac,bodyforce2,Max_Force,bondlist,BFmultiplier,Critical_ts_steel_elastic,DisplacedCoordinates,UndeformedLength)
% Calculate bond forces, make use of logical indexing

nodei=bondlist(:,1);
nodej=bondlist(:,2);

% DisplacedCoordinates=gpuArray(DisplacedCoordinates);
Xdeforme=DisplacedCoordinates(nodej,:)-DisplacedCoordinates(nodei,:);
Xdeformed(:,1)=Xdeforme(:,1);
Ydeformed(:,1)=Xdeforme(:,2);
Zdeformed(:,1)=Xdeforme(:,3);

% DeformedLength(:,1)=sqrt(Xdeforme(:,1).^2+Xdeforme(:,2).^2+Xdeforme(:,3).^2);
% Stretch=(DeformedLength-UndeformedLength)./DeformedLength;  %strain'

[DeformedLength,Stretch]=arrayfun(@gpu_subfunc2,Xdeformed,Ydeformed,Zdeformed,UndeformedLength);
Stretch=gpuArray(single(Stretch));

fail(fail==1 & BondType==0 & Stretch>Critical_ts_conc)=0;     % Deactivate bond if stretch exceeds critical stretch   Failed = 0 
fail(fail==1 & BondType==1 & Stretch>3*Critical_ts_conc)=0;   % EMU user manual recommends that the critical stretch and bond force are multiplied by a factor of 3 for concrete to steel bonds 
fail(fail==1 & BondType==2 & Stretch>Critical_ts_steel)=0;   
% fail1=gpuArray(fail);
% Bond remains active = 1
% Calculate X,Y,Z component of bond force


[BforceX,BforceY,BforceZ]=arrayfun( @gpu_subfunc1, ...
BFmultiplier,fail,const,Stretch,Volume,fac,Xdeformed,Ydeformed,Zdeformed,DeformedLength,Critical_ts_steel_elastic,Critical_ts_steel,BondType);
%check steel perfect plastic
% plastic=gpuArray(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% ultimate=gpuArray(fail==1 & BondType==2 & Stretch>0.01);
% 
% if isempty(single(plastic(plastic==true)))==0
%     ans=find(single(plastic)==1);
%     BforceX(ans)=BforceX(ans)*Critical_ts_steel_elastic./Stretch(ans);
%     BforceY(ans)=BforceY(ans)*Critical_ts_steel_elastic./Stretch(ans);
%     BforceZ(ans)=BforceZ(ans)*Critical_ts_steel_elastic./Stretch(ans);
% end
% 
% if isempty(single(ultimate(ultimate==true)))==0
%     ans=find(single(ultimate)==1);
%     BforceX(ans)=BforceX(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
%     BforceY(ans)=BforceY(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
%     BforceZ(ans)=BforceZ(ans)*2.5.*(Stretch(ans)-0.01)./(Critical_ts_steel-0.01);
% end

% [BforceX,BforceY,BforceZ] = arrayfun( @gpu_subfunc3,fail,BondType,Stretch,Critical_ts_steel_elastic,Critical_ts_steel,BforceX,BforceY,BforceZ);
% %
% % BforceX(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceX(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% BforceY(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceY(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% BforceZ(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)=BforceZ(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic)*Critical_ts_steel_elastic./Stretch(fail==1 & BondType==2 & Stretch>Critical_ts_steel_elastic);
% 
% % DeformedLength will be 0 intially and divison leads to 'Not-a-Number'. If Bforce = 'NaN', set to 0 
% BforceX(isnan(BforceX))=0; 
% BforceY(isnan(BforceY))=0;
% BforceZ(isnan(BforceZ))=0;

% 
% [BforceX1,BforceY1,BforceZ1]=arrayfun(@gpu_subfunc3,BforceX,BforceY,BforceZ,fail,BondType,Stretch,Critical_ts_steel_elastic);
% % BforceX1(isnan(BforceX1))=0; 
% % BforceY1(isnan(BforceX1))=0;
% BforceZ1(isnan(BforceX1))=0;

bondlist=gpuArray(bondlist);
BforceX2=[accumarray(bondlist(:,1),BforceX);0];
BforceY2=[accumarray(bondlist(:,1),BforceY);0];
BforceZ2=[accumarray(bondlist(:,1),BforceZ);0];
BforceX3=accumarray(bondlist(:,2),BforceX);
BforceY3=accumarray(bondlist(:,2),BforceY);
BforceZ3=accumarray(bondlist(:,2),BforceZ);
% Bond force is positive on Node i
                                   
% Add body force
% Nforce=gpuArray(Nforce);
Nforce(:,:)=Nforce+[BforceX2 BforceY2 BforceZ2]-[BforceX3 BforceY3 BforceZ3]+(bodyforce2(:,:)*Max_Force);
LOAD=sum((bodyforce2(:,:)*Max_Force));
% Nforce=gather(Nforce);
end
