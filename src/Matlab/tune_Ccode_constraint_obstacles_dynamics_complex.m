 
state =ones(8,1); 


state =[15.8828, 0, 0, 0, 0, 18.1959, 0, 0];

t = 0;

traj_ob.no  = 1; 
traj_ob.pos  = [90.11; 0.342];
traj_ob.vel= [0; 0]; 
traj_ob.acc = [0; 0]; 
traj_ob.rad = 2.1;

test_resl = constraint_obstacles_dynamics_complex(state,t,traj_ob); 


% xp_dot = u(1);  %longitudinal speed
% yp_dot = u(2);  %lateral speed
% psi_dot = u(3); 
% epsi = u(4);
% ey= u(5);  %lateral position
% s = u(6);  %logitudinal position 
% 
% %20181011 add: 
% steer = u(7);   %steering angle 
% acc = u(8);  %longitudinal acc 
% 
%  
% %static obstacles: 
% % global no_ob;
% global dead;  %once there is vilation of constraints, trigger dead,
% % global pos_ob_array_pre radius_pre;
% no_ob_input = traj_ob.no;
% pos_ob_array_pre_input = zeros(2,no_ob_input);
% vel_ob_array_pre = zeros(2,no_ob_input);
% acc_ob_array = zeros(2,no_ob_input);
% 
% pos_ob_array_pre_input = traj_ob.pos;
% vel_ob_array_pre = traj_ob.vel;
% acc_ob_array = traj_ob.acc;
% radius_pre = traj_ob.rad; 