function [] = calculateshear(nDIVX,DX,COORDINATES,stresstensor)
% Calculate the shear force along the length (x-axis) of the member

datageometry
nodesPlaneYZ=zeros(nDIVX,1);
shearForce=zeros(nDIVX,1);

for i=1:nDIVX % Loop the member length
    
    % Identify and flag all nodes within cross-sectional plane Y-Z
    nodesPlaneYZ(i,1)=i*DX;
    crossSectionFlag=(COORDINATES(:,1)==nodesPlaneYZ(i,1))==1; 
    coordCrossSection=COORDINATES(:,:);
    stressCrossSection=stresstensor(:,:,:);
    logicCondition1 = crossSectionFlag==0;          % Delete node if it is not located in cross-section (flag==0)
    coordCrossSection(logicCondition1,:)=[];
    stressCrossSection(logicCondition1,:)=[];       % TODO: I believe that this line re-organises the 3x3 stress matrix into a 1x9 row. Check if this is true.
   
    shearForce(i,1)=(sum(stressCrossSection(:,3))/size(coordCrossSection,1))*(WIDTH*HEIGHT); % Shear Force = average shear stress / cross-sectional area
    
end

plot(nodesPlaneYZ,-shearForce);
str1=sprintf('Shear Force Diagram, average shear force = %.0fN', mean(-shearForce));
title(str1)
xlabel('x axis (m)')
ylabel('Shear Force V (N)')

end

