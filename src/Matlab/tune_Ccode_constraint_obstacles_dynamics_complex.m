 
% state =ones(8,1); 


state =[15.8828, 0.6, 0.1, 2, 2.3, 18.1959, 0, 0];

t = 0;

%test data, the same with the data in C++:
traj_ob.no  = 5; 
traj_ob.pos(:,1)  = [108.233; -1.36608];
traj_ob.pos(:,2)  = [123.244; 0.998566];
traj_ob.pos(:,3)  = [142.236; -0.472505];
traj_ob.pos(:,4)  = [157.027; -1.60909];
traj_ob.pos(:,5)  = [169.361; -0.807626];
traj_ob.vel= zeros(2,traj_ob.no);
traj_ob.acc = zeros(2,traj_ob.no);
traj_ob.rad = [2.00564; 2.62428; 2.23689; 2.28824; 3.25264];

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


% 
% L_f_h_ang_part1 = - ((psi_dot - psi_dot_com)*(xp_dot^2 + yp_dot^2 - vel_ob_x*xp_dot*cos(epsi) - ...
%     vel_ob_y*yp_dot*cos(epsi) - vel_ob_y*xp_dot*sin(epsi) + vel_ob_x*yp_dot*sin(epsi))*(ey*vel_ob_x ...
%     + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)) - ((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)) - ((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)) - ((vel_ob_x*cos(epsi) - xp_dot + vel_ob_y*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2));
        