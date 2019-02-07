
function [Bond_damage]=Damage(perc_progress,tt0,NumFamMembVector,fail,Totalnodes,Totalbonds,coordinates,disp,bondlist,CrossSectionFlag)

% Damage - Calculate the damage (percentage of bonds broken) for every node

% Initialise
Bond_damage=gpuArray.zeros(Totalnodes,1,'single');
UnbrokenBonds=gpuArray.zeros(Totalnodes,1,'single');

% Calculate the number of unbroken bonds attached to every node
UnbrokenBonds=[accumarray(bondlist(:,1),fail);0]+UnbrokenBonds;
UnbrokenBonds=accumarray(bondlist(:,2),fail)+UnbrokenBonds;

Bond_damage(:,1)=1-(UnbrokenBonds./NumFamMembVector(:,1)); % Calculate percentage of broken bonds attached to node i (1 would indicate that all bonds have broken)

% Plot data
% DSF=1; % Displacement scale factor
% pointsize=1;
% figure(1);
% %subplot(1,2,1)
% scatter3(coordinates(:,1)+(disp(:,1)*DSF),coordinates(:,2)+(disp(:,2)*DSF),coordinates(:,3)+(disp(:,3)*DSF),pointsize,Bond_damage(:,1));
% axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% %view([0,-90,0]) % View in 2D
% view(30,30);    % View in 3D
% grid off;
% colormap jet ;
% %colorbar
% caxis([0 1]);
% h = colorbar;
% ylabel(h, 'Damage')
% title(num2str(perc_progress));
% drawnow;
% 
% filename = 'damage3d.gif';
% frame = getframe(1);
% im = frame2im(frame);
% [A,map] = rgb2ind(im,256);
% if tt0 == 1
%     imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);    
% else    
%     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);    
% end
% saveas(1,['damage_' num2str(perc_progress) '.fig'],'fig');
%  
% % % Plot cross section of data
% % % Seperate scatter data into sub sets
% 
% 
% CoordCrossSection=coordinates(:,:);
% DispCrossSection=disp(:,:);
% BondDamageCrossSection=Bond_damage(:,:);
% LogicCondition1 = CrossSectionFlag==0; % Delete node if it is not located in cross-section
% CoordCrossSection(LogicCondition1,:)=[];
% DispCrossSection(LogicCondition1,:)=[];
% BondDamageCrossSection(LogicCondition1,:)=[];
% 
% pointsize=1;
% figure(2);
% scatter(CoordCrossSection(:,1)+(DispCrossSection(:,1,1)*DSF),CoordCrossSection(:,3)+(DispCrossSection(:,3,1)*DSF),pointsize,BondDamageCrossSection(:,1))
% axis equal;
% xlabel('x');
% ylabel('z');
% colormap jet ;
% colorbar;
% caxis([0 1]);
% h = colorbar;
% ylabel(h, 'Damage');
% title(num2str(perc_progress));
% drawnow;
% % 
% % 
% filename = 'damage2d.gif';
% frame = getframe(2);
% im = frame2im(frame);
% [A,map] = rgb2ind(im,256);
% if tt0 == 1
%     imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);    
% else    
%     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);    
% end

end 