function [DeformedLength,Xdeformed,Ydeformed,Zdeformed,Stretch]=DeformedLengthfunc(Totalbonds,bondlist,UndeformedLength,DeformedLength,coordinates,disp,Xdeformed,Ydeformed,Zdeformed)

% Calculate the deformed length of a bond using a nested for loop
DisplacedCoordinates=coordinates+disp;

% for i=1:Totalbonds                           
%     
%    nodei=bondlist(i,1);
%    nodej=bondlist(i,2);
% 
%     % Length of deformed bond
%     DeformedLength(i)=(DisplacedCoordinates(nodej,1)-DisplacedCoordinates(nodei,1))^2+(DisplacedCoordinates(nodej,2)-DisplacedCoordinates(nodei,2))^2+(DisplacedCoordinates(nodej,3)-DisplacedCoordinates(nodei,3))^2;
%     Xdeformed(i)=(DisplacedCoordinates(nodej,1)-DisplacedCoordinates(nodei,1)); % X-component of deformed bond
%     Ydeformed(i)=(DisplacedCoordinates(nodej,2)-DisplacedCoordinates(nodei,2)); % Y-component of deformed bond
%     Zdeformed(i)=(DisplacedCoordinates(nodej,3)-DisplacedCoordinates(nodei,3)); % Z-component of deformed bond
% 
% end
% DeformedLength=sqrt(DeformedLength);

nodei=bondlist(:,1);
nodej=bondlist(:,2);
Xdeforme=DisplacedCoordinates(nodej,:)-DisplacedCoordinates(nodei,:);
DeformedLength(:,1)=sqrt(Xdeforme(:,1).^2+Xdeforme(:,2).^2+Xdeforme(:,3).^2);
Xdeformed(:,1)=Xdeforme(:,1);
Ydeformed(:,1)=Xdeforme(:,2);
Zdeformed(:,1)=Xdeforme(:,3);
Stretch=(DeformedLength-UndeformedLength)./DeformedLength;  %strain'

end

