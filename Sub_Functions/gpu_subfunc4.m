function [BforceX,BforceY,BforceZ] = gpu_subfunc4(fail,BondType,Stretch,Critical_ts_steel_elastic,Critical_ts_steel,BforceX,BforceY,BforceZ)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

while fail==1 && BondType==2 && Stretch>Critical_ts_steel_elastic
    BforceX=BforceX*Critical_ts_steel_elastic/Stretch;
    BforceY=BforceY*Critical_ts_steel_elastic/Stretch;
    BforceZ=BforceZ*Critical_ts_steel_elastic/Stretch;
end

while fail==1 && BondType==2 && Stretch>0.01
    BforceX=BforceX*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
    BforceY=BforceY*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
    BforceZ=BforceZ*2.5*(Stretch-0.01)/(Critical_ts_steel-0.01);
end



% while fail==false
%     BforceX=single(0);
%     BforceY=single(0);
%     BforceZ=single(0);
% end

end

