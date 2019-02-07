% Plot displacement vs time

function ReactionForceVsDisplacement(ReactionForce,Displacement_node,LOAD_step,dx,failnum)


% figure(6);
% plot(-Displacement_node,-ReactionForce)
% xlabel('Displacement (m)')
% ylabel('ReactionForce') 

figure(6);
set(gcf, 'Position', 0.9*get(0, 'Screensize'));

subplot(3,1,1);
plot([ReactionForce(:,1) LOAD_step(:,1)]*dx*1e-6/6.895/2);
xlabel('Step');
ylabel('X-axis (ksi)') ;
legend({'Reaction','Load'},'Location','southeast');

% subplot(3,1,2);
% plot([ReactionForce(:,2)*dx^3 LOAD_step(:,2)*dx^3]*dx*1e-6/6.895/2);
% % xlabel('Step')
%  ylabel('Norm Y-axis (ksi)') ;
% 
% subplot(3,1,3);
% plot([ReactionForce(:,3)*dx^3 LOAD_step(:,3)*dx^3]*dx*1e-6/6.895/2);
% % xlabel('Step')
%  ylabel('Norm Z-axis (ksi)') ;

subplot(3,1,2);
plot(ReactionForce(:,2:3)*dx*1e-6/6.895/2);
% xlabel('Step')
ylabel('Norm Y/Z-axis (ksi)') ;
legend({'Y-axis','Z-axis'},'Location','southeast');

subplot(3,1,3);
plot(failnum);
% xlabel('Step')
 ylabel('Num of Fail') ;

drawnow;
end 