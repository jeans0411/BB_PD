function [DeformedLength,Stretch] = gpu_subfunc2(Xdeformed,Ydeformed,Zdeformed,UndeformedLength)

DeformedLength=sqrt(Xdeformed.^2+Ydeformed.^2+Zdeformed.^2);
Stretch=(DeformedLength-UndeformedLength)/DeformedLength;  %strain'

end

