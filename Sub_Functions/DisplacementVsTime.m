% Plot displacement vs time

function DisplacementVsTime(nt,Displacement_node,LOAD_step,dx)

% Plot y displacement of bottom node at the free end - node 200
%Displacement_node=-squeeze(Displacement(200,2,:));
% 
Time(:,1)=1:nt;

figure;
plot(Displacement_node)
xlabel('Time step') 
ylabel('Displacement (m)')

figure;
plot(Displacement_node*1e6,LOAD_step*dx*1e-6/6.895/2,'*')
xlabel('Strain (us)') 
ylabel('Rebar Stress (ksi)')

end 