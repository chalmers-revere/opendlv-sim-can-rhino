%% load the data, modified on 20190919

% load  data_nominal_states.txt;
% load data_msg_nom_u.txt;
% load  data_model_state.txt;
% load data_msg_actual_u.txt;
% load data_ref_traj.txt;

clear all;

filename=dir('data_nominal_states_*.txt');
data_nominal_states = load(filename.name); 

filename=dir('data_msg_nom_u_*.txt');
data_msg_nom_u = load(filename.name); 

filename=dir('data_model_state_*.txt');
data_model_state = load(filename.name); 

filename=dir('data_msg_actual_u_*.txt');
data_msg_actual_u = load(filename.name); 

filename=dir('data_ref_traj_*.txt');
data_ref_traj = load(filename.name); 

u_ctrl = data_msg_nom_u(:,2:3); %nominal 
u_actual = data_msg_actual_u(:,2:3);  %actual 

y1_nom = data_nominal_states(:,2:9);
T_sampl_nom = 1/150; 
t = T_sampl_nom:T_sampl_nom:T_sampl_nom*size(y1_nom,1);

y1_actual = data_model_state(:,2:9);
t_actual = T_sampl_nom:T_sampl_nom:T_sampl_nom*size(data_model_state,1);
t_uactual = T_sampl_nom:T_sampl_nom:T_sampl_nom*size(u_actual,1);

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
filename=dir('data_traj_ob_*.txt');
data_traj_ob = load(filename.name); 
 
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


if length(delta)< length(dot_y)
    delta = [delta; delta(end)* (ones(length(dot_y) - length(delta), 1))];
elseif length(delta)> length(dot_y)
%     dot_y = [dot_y; dot_y(end)* (ones(length(delta) - length(dot_y), 1))];
%     dot_x = [dot_x; dot_x(end)* (ones(length(delta) - length(dot_x), 1))];
%     dot_psi = [dot_psi; dot_psi(end)* (ones(length(delta) - length(dot_psi), 1))];
    delta = delta(1:length(dot_y));
end
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
plot(t, data_ref_traj(:,3), '--', t,P_nom(:,1),'-.',t_actual,P_sens(:,1)),grid;
ylabel('X(m)');
title('POSITION:SENSED VS COMMAND');
 

subplot(4,1,2);
plot(t, data_ref_traj(:,2), '--',t,P_nom(:,2),'-.',t_actual,P_sens(:,2)),grid;
ylabel('Y(m)');
 xlabel('time(s)');
 
 subplot(4,1,3);
plot(t, data_ref_traj(:,6), '--',  t, y1_nom(:,1), '-.',  t_actual,y1_actual(:,1) ),grid;
legend('reference','nominal', 'actual');
ylabel('v(m/s)');
 xlabel('time(s)');
 
 subplot(4,1,4);
plot( t, y1_nom(:,4)/pi*180, '-.',  t_actual,y1_actual(:,4)/pi*180 ),grid;
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
  plot(t_uactual, u_actual(:,1));  title('actual control'); 
ylabel('  a_x(m/s^2)');
 subplot(2,1,2);
  plot(t_uactual, u_actual(:,2)); 
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
 
 load data_safety_certificate.txt;
 
 An_sidepos = data_safety_certificate(:,1:2);
 bn_sidepos= data_safety_certificate(:,3); 
 
 %notice the order: 
 test_cons = An_sidepos(:, 1).*u_ctrl(:,2) + An_sidepos(:, 2).*u_ctrl(:,1) - bn_sidepos; 
 
 figure(8);
   subplot(5,1,1);
  plot(t, data_safety_certificate(:,1));  grid;  title('A1, side'); 
    subplot(5,1,2);
  plot(t, data_safety_certificate(:,2));  grid;  title('A2, side'); 
    subplot(5,1,3);
  plot(t, data_safety_certificate(:,3));  grid;  title('b, side'); 
  
      subplot(5,1,4);
  plot(t, test_cons);  grid;  title('test_cons, pos side'); 
  

  
  subplot(5,1,5);
  plot(t, data_safety_certificate(:,4));  grid;  title('cbf, side'); 
  
figure(9);
  plot(t, data_safety_certificate(:, 9));  grid;  title('min_value, qp'); 
  
% % ylabel('  a_x(m/s^2)');
%  subplot(2,1,2);
%   plot(t, data_safety_certificate(:,8));  grid;  
% % ylabel('  steer(rad)');
 xlabel('time(s)');
 
 