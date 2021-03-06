function [nodefamily,nfpointer,UndeformedLength,NumFamMembVector,NumFamMembVector2,bondlist,bondlist2,Totalbonds]=BuildHorizons(Totalnodes,coordinates,delta)

% Determine the nodes inside the horizon of each material point, build bond lists, and determine undeformed length of every bond
%bondlist is the core.
%% KD-tree
[NF]=rangesearch(coordinates,coordinates,delta); % rangesearch is the KDTreeSearcher function for distance search.


%% Create node family data structure
counter=1;
for i=1:Totalnodes
    
    CurrentNodeFamilyTemp=NF{i}; % First element indicates the origin node
    CurrentNodeFamilyTemp(1)=[]; % Remove first element (origin node)
    %CurrentNodeFamily=sort(CurrentNodeFamilyTemp); % Sort nodes in ascending order
    CurrentNodeFamily=CurrentNodeFamilyTemp; 
    
    NumFamMemb=size(CurrentNodeFamily,2);
    NumFamMembVector(i,1)=NumFamMemb;    
    
    % Create list of pointers to node family members
    if i==1
        nfpointer(i,1)=1;
    else
        nfpointer(i,1)=nfpointer(i-1)+NumFamMembVector(i-1);
    end
      
    for j=1:NumFamMemb
        nodefamily(counter,1)=CurrentNodeFamily(1,j);
        counter=counter+1;
    end
    
end

%% Create bond list

bondlist=zeros(size(nodefamily,1)/2,2);
counter1=0;
NumFamMembVector2=NumFamMembVector;
for i=1:Totalnodes
    % All nodes within Node 'i' sphere of influence
    for j=1:NumFamMembVector(i)
        % Consider bond between Node 'i' and Node 'cnode'
        cnode=nodefamily(nfpointer(i)+(j-1),1);
        if cnode>i % if i is greater than cnode, the corresponding bond has already been added to the list. The following link explains other methods to prevent the double checking of indices https://stackoverflow.com/questions/31961009/java-how-to-stop-nested-loops-from-checking-same-indices-twice
            counter1=counter1+1;
            bondlist(counter1,:)=[i cnode];
        else
            NumFamMembVector2(i)=NumFamMembVector2(i)-1;
        end
    end
end
Totalbonds=counter1;
% Totalbonds=size(bondlist,1);

%% Calculate undeformed length of every bond

UndeformedLength=zeros(Totalbonds,1);

for i=1:size(bondlist,1)
    nodei=bondlist(i,1);
    nodej=bondlist(i,2);
    UndeformedLength(i)=(coordinates(nodei,1)-coordinates(nodej,1))^2+(coordinates(nodei,2)-coordinates(nodej,2))^2+(coordinates(nodei,3)-coordinates(nodej,3))^2;
end

UndeformedLength=sqrt(UndeformedLength);
tic;
bondlist2=nan(Totalnodes,max(NumFamMembVector2));
i=1;j1=1;
while i<=Totalnodes
    bondlist2(i,1:NumFamMembVector2(i))=bondlist(j1:j1+NumFamMembVector2(i)-1,2)';
    j1=NumFamMembVector2(i)+j1;
    i=i+1;
end
toc
end