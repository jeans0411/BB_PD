function [BforceX,BforceY,BforceZ] = gpu_subfunc1(BFmultiplier,fail,const,Stretch,Volume,fac,Xdeformed,Ydeformed,Zdeformed,DeformedLength,Critical_ts_steel_elastic,Critical_ts_steel,BondType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
BforceX=BFmultiplier*fail*const*Stretch*Volume*fac*(Xdeformed/DeformedLength);
BforceY=BFmultiplier*fail*const*Stretch*Volume*fac*(Ydeformed/DeformedLength);
BforceZ=BFmultiplier*fail*const*Stretch*Volume*fac*(Zdeformed/DeformedLength);

if fail==1 && BondType==2 && Stretch>Critical_ts_steel_elastic
    BforceX=BforceX*Critical_ts_steel_elastic/Stretch;
    BforceY=BforceY*Critical_ts_steel_elastic/Stretch;
    BforceZ=BforceZ*Critical_ts_steel_elastic/Stretch;
end

if fail==1 && BondType==2 && Stretch>0.01
    BforceX=BforceX*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
    BforceY=BforceY*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
    BforceZ=BforceZ*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
end
end

