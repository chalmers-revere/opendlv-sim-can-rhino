%load the data 
load  data_nominal_states.txt;
load data_msg_nom_u.txt;
load  data_model_state.txt;
load data_msg_actual_u.txt;


u_ctrl = data_msg_nom_u(:,2:3); %nominal 
u_actual = data_msg_actual_u(:,2:3);  %actual 

y1_nom = data_nominal_states(:,2:9);
T_sampl_nom = 0.01; 
t = T_sampl_nom:T_sampl_nom:T_sampl_nom*size(y1_nom,1);

y1_actual = data_model_state(:,2:9);
t_actual = T_sampl_nom:T_sampl_nom:T_sampl_nom*size(data_model_state,1);


figure(1);
subplot(4,1,1);
plot(t , y1_nom(:,1), '-.', t_actual,  y1_actual(:,1));
ylabel('v_x');

subplot(4,1,2);
plot(t , y1_nom(:,2), '-.', t_actual,  y1_actual(:,2));
ylabel('v_y');

subplot(4,1,3);
plot(t , y1_nom(:,5), '-.', t_actual,  y1_actual(:,5));
ylabel(  'e_y'); 
xlabel('Time');

subplot(4,1,4);
plot(t , y1_nom(:,6), '-.', t_actual,  y1_actual(:,6)); grid; 
ylabel(  's'); 
legend('nominal', 'actual');
xlabel('Time');  




close all;


% load sim_data.mat;


% t = t1;

P_nom =   y1_nom(:, [6,5]); 
% y1_actual = y1_nom;
P_sens = y1_actual(:, [6,5]); 
P_cente = P_sens;

%%0815, record the reference trajectory: 
% reference = u_ctrl(3:end,:);   %9-by-n, traj, dot traj, ddot traj,  here traj = [phi, ey, s];
% P_cente=reference([3,2],:)';
%0815, traj_ob_seris: 2*n-by-no_ob

 load data_traj_ob.txt;
%  traj_ob_seris = data_traj_ob; 
 
traj_ob_seris = data_traj_ob;
%0815, traj_ob_seris: 2*n-by-no_ob
no_ob = size(traj_ob_seris, 2);
traj_ob_plot = zeros(2, size(traj_ob_seris,1)/3, no_ob);
for i_ob =1:no_ob
    for i_time =1:size(traj_ob_seris,1)/3
        traj_ob_plot(:,i_time,i_ob) = traj_ob_seris((i_time-1)*3+1:(i_time-1)*3+2,i_ob);
    end
end

% 
% traj_ob_plot = zeros(2, length(t1), no_ob);
% for i_ob =1:no_ob
%     for i_time =1:length(t1)
%         traj_ob_plot(:,i_time,i_ob) = traj_ob_seris((i_time-1)*2+1:(i_time-1)*2+2,i_ob);
%     end
% end


%%road side:
road_side_x = min(P_sens(:,1)):0.1:max(P_sens(:,1));
road_side_y1 = 3.7*ones(1, length(road_side_x));
road_side_y2 = -3.7*ones(1, length(road_side_x));

delta = u_actual(:,2); 
% delta = 0; 
dot_y = y1_actual(:,2);
dot_x = y1_actual(:,1);
dot_psi = y1_actual(:,3); 

a = 1.68; 
b = 1.715; 

 
alpha_f= (dot_y+a*dot_psi)./dot_x - delta; 
alpha_r = (dot_y -b*dot_psi)./dot_x; 
alpha_f(dot_x<=1e-1) = 0;
alpha_r(dot_x<=1e-1) = 0;


figure(100); 
subplot(2,1,1);
plot(t_actual,alpha_f ),grid;
ylabel('\alpha_f');
title('slip ratios');
% legend('P_c_o_m','P_s_e_n');

subplot(2,1,2);
plot(t_actual,alpha_r ),grid;
ylabel('\alpha_r');
 xlabel('time(s)');


len = length(t_actual); 
h_plot = zeros(len,1);

for i=1:len
    h_plot(i) = (P_sens(i,1) - 80)^2 + (P_sens(i,2) - 0.5)^2 -1;
end
 

figure(1); 
%the tube along nominal trajectory: 
dm = 1.414; %maximum disturbance 
k1 = 3; 
r_tube = dm/k1; 
for ii = 1:10:length(t)     
    circle_tube(r_tube,P_nom(ii,1), P_nom(ii,2)); hold on;
end


% figure(1);   
plot(P_cente(:,1),P_cente(:,2),'-.', P_nom(:,1),P_nom(:,2), '--', P_sens(:,1),P_sens(:,2)),grid; hold on; legend('c','p_d', 'p');
plot(road_side_x,road_side_y1,'-.k',road_side_x,road_side_y2,'-.k'); hold on;
% axis equal;
xlabel('X(m)');ylabel('Y(m)');
title('POSITION TRACKING:SENSED VS COMMAND');
axis equal; 

% global pos_ob_array_pre radius_pre;
%pos_ob_array_pre(2, :) = - pos_ob_array_pre(2,:);  %only for plot the mpc data with mistake 
for i=1:no_ob
    plot(traj_ob_plot(1,1,i), traj_ob_plot(2, 1,i), '*r'); hold on; 
    Ds = traj_ob_seris(3, i);
    circle(Ds,traj_ob_plot(1,1,i), traj_ob_plot(2, 1,i)); hold on; 
end

%the trajectory of the obstacles: 
for i_ob =1:no_ob    
    plot(traj_ob_plot(1,:,i_ob), traj_ob_plot(2,:,i_ob), 'LineWidth', 4); hold on; 
end
 
 
figure(2);   
subplot(4,1,1);
plot(t,P_nom(:,1),'-.',t_actual,P_sens(:,1)),grid;
ylabel('X(m)');
title('POSITION:SENSED VS COMMAND');
 

subplot(4,1,2);
plot(t,P_nom(:,2),'-.',t_actual,P_sens(:,2)),grid;
ylabel('Y(m)');
 xlabel('time(s)');
 
 subplot(4,1,3);
plot(  t, y1_nom(:,1), '--',  t_actual,y1_actual(:,1) ),grid;
legend('reference','nominal', 'actual');
ylabel('v(m/s)');
 xlabel('time(s)');
 
 subplot(4,1,4);
plot( t, y1_nom(:,4)/pi*180, '--',  t_actual,y1_actual(:,4)/pi*180 ),grid;
ylabel('\psi(degree)');
 xlabel('time(s)');
 
 figure(20); 
  subplot(2,1,1);
plot(  t, y1_nom(:,2), '--',  t_actual,y1_actual(:,2) ),grid;
legend( 'nominal', 'actual');
ylabel('v_y(m/s)');
 xlabel('time(s)');
 
 subplot(2,1,2);
plot(  t, y1_nom(:,3), '--',  t_actual,y1_actual(:,3) ),grid;
ylabel('\dot \psi(rad/s)');
 xlabel('time(s)');
 
 
%  figure(3);
%  plot(t, h_plot); 
% ylabel('distance');
%  xlabel('time(s)');
 
 
 
 figure(4); 
 subplot(2,1,1);
  plot(t, u_ctrl(:,1));  title('nominal control'); 
ylabel(' a_x(m/s^2)');
 subplot(2,1,2);
  plot(t, u_ctrl(:,2)); 
ylabel(' steer(rad)');
 xlabel('time(s)');
 
 
 figure(5); 
 subplot(2,1,1);
  plot(t_actual, u_actual(:,1));  title('actual control'); 
ylabel('  a_x(m/s^2)');
 subplot(2,1,2);
  plot(t_actual, u_actual(:,2)); 
ylabel('  steer(rad)');
 xlabel('time(s)');
 
 
  figure(6); 
 subplot(2,1,1);
  plot(t, y1_nom(:,8));  title('nominal control, filtered'); 
ylabel('  a_x(m/s^2)');
 subplot(2,1,2);
  plot(t, y1_nom(:,7)); 
ylabel('  steer(rad)');
 xlabel('time(s)');
 
 
   figure(7); 
 subplot(2,1,1);
  plot(t_actual, y1_actual(:,8));  title('actual control, filtered'); 
ylabel('  a_x(m/s^2)');
 subplot(2,1,2);
  plot(t_actual, y1_actual(:,7)); 
ylabel('  steer(rad)');
 xlabel('time(s)');
 
 