clear all; 
close all;

% QUESTIONS:
% is delta_u(t-delta_t) correct for first itteration of for loop (i.e. its set to zero)?
% do we impliment the 4mm check ourselves or not for SimulatorScriptOfficial?
% can only get within 4cm

m1 = 0.375;
m2 = 0.375;
l1 = 0.15;
l2 = 0.15;
c1 = 2;
c2 = 2;
g = 3.7;

lin_space_points = 5;
delta_T = 0.001;
t_max = 2;
tspan = 0:delta_T:t_max;

% initializing all variables for each via point
A = sym(zeros(4,4,lin_space_points*4-4));
B = sym(zeros(4,2,lin_space_points*4-4));
C = sym(zeros(2,4,lin_space_points*4-4));
D = sym(zeros(2,2,lin_space_points*4-4));

ev_A = zeros(4,1,lin_space_points*4-4);
ev_A = complex(ev_A);

K = zeros(2,4,lin_space_points*4-4);
K = complex(K);

F = zeros(4,2,lin_space_points*4-4);
F = complex(F);

x_op = zeros(4,1,lin_space_points*4-4);
tau1_op = zeros(lin_space_points*4-4,1);
tau2_op = zeros(lin_space_points*4-4,1);
u_op = zeros(2,1,lin_space_points*4-4);

delta_x_hat = zeros(4,1,length(tspan));
u_input = zeros(2,1,length(tspan));
u_op_each_itteration = zeros(2,1,length(tspan));

q1_visualize = zeros(size(tspan,2),1);
q2_visualize = zeros(size(tspan,2),1);
q1d_visualize = zeros(size(tspan,2),1);
q2d_visualize = zeros(size(tspan,2),1);

q1_via_visualize = zeros(size(tspan,2),1);
q2_via_visualize = zeros(size(tspan,2),1);

Q_lqr = zeros(4,4,lin_space_points*4-4);
R_lqr = zeros(2,2,lin_space_points*4-4);
Q_kf = zeros(4,4,lin_space_points*4-4);
R_kf = zeros(2,2,lin_space_points*4-4);

mean = 0;
std_dev = 1/3*(pi/180);
variance = std_dev^2;

%% IK Model

x_coord = [0.1 0.2 0.2 0.1];
y_coord = [0.2 0.2 0.1 0.1];

q1_angle = zeros(1,4);
q2_angle = zeros(1,4);

for i=1:4
    q1_angle(i) = atan2(y_coord(i),x_coord(i))-acos((x_coord(i)^2 + y_coord(i)^2 + l1^2 - l2^2)/(2*l1*sqrt(x_coord(i)^2+y_coord(i)^2)));
    q2_angle(i) = acos((x_coord(i)^2+y_coord(i)^2-l1^2-l2^2)/(2*l1*l2));
end

%% FK to verify

for i=1:4
    x_FK(i) = l1*cos(q1_angle(i))+l2*cos(q1_angle(i)+q2_angle(i));
    y_FK(i) = l1*sin(q1_angle(i))+l2*sin(q1_angle(i)+q2_angle(i));
end

%% Total angular displacement

TAD = abs(q1_angle(4)-q1_angle(3))+abs(q1_angle(3)-q1_angle(2))+abs(q1_angle(2)-q1_angle(1))+abs(q1_angle(1)-q1_angle(4))+ abs(q2_angle(4)-q2_angle(3))+abs(q2_angle(3)-q2_angle(2))+abs(q2_angle(2)-q2_angle(1))+abs(q2_angle(1)-q2_angle(4));

%% EOAT Position plot

q1_via_point(1,:) = linspace(q1_angle(1),q1_angle(2),lin_space_points);
q1_via_point(2,:) = linspace(q1_angle(2),q1_angle(3),lin_space_points);
q1_via_point(3,:) = linspace(q1_angle(3),q1_angle(4),lin_space_points);
q1_via_point(4,:) = linspace(q1_angle(4),q1_angle(1),lin_space_points);

q2_via_point(1,:) = linspace(q2_angle(1),q2_angle(2),lin_space_points);
q2_via_point(2,:) = linspace(q2_angle(2),q2_angle(3),lin_space_points);
q2_via_point(3,:) = linspace(q2_angle(3),q2_angle(4),lin_space_points);
q2_via_point(4,:) = linspace(q2_angle(4),q2_angle(1),lin_space_points);

x_plot_via_point = [l1*cos(q1_via_point(1,:))+l2*cos(q1_via_point(1,:)+q2_via_point(1,:)),l1*cos(q1_via_point(2,:))+l2*cos(q1_via_point(2,:)+q2_via_point(2,:)),l1*cos(q1_via_point(3,:))+l2*cos(q1_via_point(3,:)+q2_via_point(3,:)),l1*cos(q1_via_point(4,:))+l2*cos(q1_via_point(4,:)+q2_via_point(4,:))];
y_plot_via_point = [l1*sin(q1_via_point(1,:))+l2*sin(q1_via_point(1,:)+q2_via_point(1,:)),l1*sin(q1_via_point(2,:))+l2*sin(q1_via_point(2,:)+q2_via_point(2,:)),l1*sin(q1_via_point(3,:))+l2*sin(q1_via_point(3,:)+q2_via_point(3,:)),l1*sin(q1_via_point(4,:))+l2*sin(q1_via_point(4,:)+q2_via_point(4,:))];

% converting to single row vectors
q1_via_point = [q1_via_point(1,(2:end)) q1_via_point(2,(2:end)) q1_via_point(3,(2:end)) q1_via_point(4,(2:end))];
q2_via_point = [q2_via_point(1,(2:end)) q2_via_point(2,(2:end)) q2_via_point(3,(2:end)) q2_via_point(4,(2:end))];

%% SS

syms q1 q2 q1_d q2_d q1_dd q2_dd tau1 tau2
eqns = [tau1 == ((m1*l1^2)/3 + (m2*l2^2)/12 + m2*(l1^2+(l2^2)/4+l1*l2*cos(q2)))*q1_dd + ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q2_dd - m2*l1*l2*sin(q2)*q1_d*q2_d - (m2*l1*l2*sin(q2))/2*q2_d^2 + ((m1*l1)/2+m2*l1)*g*cos(q1) + (m2*l2)/2*g*cos(q1+q2) + c1*q1_d , tau2 == ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q1_dd + (m2*l2^2)/3*q2_dd + (m2*l1*l2*sin(q2))/2*q1_d^2 + (m2*l2)/2*g*cos(q1+q2) + c2*q2_d];
dd_vals = solve(eqns, [q1_dd q2_dd]);
tau_vals = solve(eqns, [tau1 tau2]);

x = [q1; q1_d; q2; q2_d];
x_dot = [q1_d; dd_vals.q1_dd; q2_d; dd_vals.q2_dd];
u = [tau1; tau2];
y = [q1; q2];

i = 1;
ready_for_next_via_point = 1;
first_itteration_no_loop = 1;
first_itteration_for_loop = 1;
holding_time = 0; % time counter for oscillations to die down
dist = 0; % current distance between actual vs desired position
required_dist = 0.05; % this should be the 4mm distance

if first_itteration_no_loop == 1
    x_0 = [q1_via_point(end); 0; q2_via_point(end); 0]; % Note: 'end' is point A
    tau1_0 = double(subs(tau_vals.tau1, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(end) 0 0 q2_via_point(end) 0 0]));
    tau2_0 = double(subs(tau_vals.tau2, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(end) 0 0 q2_via_point(end) 0 0]));
    tau_0 = [tau1_0;tau2_0]; % arbitrarily chosen
    X0=x_0;
    U=tau_0;
    [tout,qout] = ode45(@(time,x)simulatorofficial(time,x,U,l1,l2,m1,m2,g,c1,c2),[0 0.001],X0);
    q=qout(end,[1,3])';
    q1_visualize(1) = q(1);
    q2_visualize(1) = q(2);
    q1d_visualize(1) = qout(end,2);
    q2d_visualize(1) = qout(end,4);
    q1_via_visualize(1) = q1_via_point(end);
    q2_via_visualize(1) = q2_via_point(end);
    first_itteration_no_loop = 0;
end

for t=0.001:0.001:t_max
    t_ind = round(t/delta_T);
    if ready_for_next_via_point == 1
        A(:,:,i) = jacobian(x_dot,x);
        B(:,:,i) = jacobian(x_dot,u);
        C(:,:,i) = jacobian(y,x);
        D(:,:,i) = jacobian(y,u);

        x_op(:,:,i) = [q1_via_point(i); 0; q2_via_point(i); 0];
        tau1_op(i) = double(subs(tau_vals.tau1, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(i) 0 0 q2_via_point(i) 0 0]));
        tau2_op(i) = double(subs(tau_vals.tau2, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(i) 0 0 q2_via_point(i) 0 0]));
        u_op(:,:,i) = [tau1_op(i); tau2_op(i)];

        % Substituing the equilibrium points and converting from syms to numbers type 
        A(:,:,i) = double(subs(A(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));       
        B(:,:,i) = double(subs(B(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
        C(:,:,i) = double(subs(C(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
        D(:,:,i) = double(subs(D(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
       
%         Q_lqr(:,:,i) = [1000 0 0 0;
%             0 10 0 0;
%             0 0 1 0;
%             0 0 0 1];
%         R_lqr(:,:,i) = [10000 0;
%             0 1];
%         [K_lqr, P_lqr, ev_lqr] = lqr(double(A(:,:,i)), double(B(:,:,i)), Q_lqr(:,:,i), R_lqr(:,:,i));
%         K(:,:,i) = K_lqr;
        
        Q_kf(:,:,i) = eye(4)*variance;
        R_kf(:,:,i) = eye(2)*variance;
        [F_transpose_kf, P_kf, ev_kf] = lqr(double(A(:,:,i)).', double(C(:,:,i)).', Q_kf(:,:,i), R_kf(:,:,i)); % Note: this returns F transpose
        F(:,:,i) = transpose(F_transpose_kf);

        K(:,:,i) = place(double(A(:,:,i)), double(B(:,:,i)), [-100 -200 -300 -400]);
%       F(:,:,i) = transpose(place(double(transpose(A(:,:,i))), double(transpose(C(:,:,i))), [-10 -21 -32 -43]));
    end
    
    u_op_each_itteration(:,:,t_ind) = u_op(:,:,i); % store value of u_op at current t_ind
    ready_for_next_via_point = 0;
    
    if first_itteration_for_loop == 1
        delta_x_hat(:,:,t_ind) = (x_0 - x_op(:,:,i)) + delta_T*((A(:,:,i)-F(:,:,i)*C(:,:,i))*(x_0 - x_op(:,:,i)) + B(:,:,i)*(U - tau_0) + F(:,:,i)*(q - [q1_via_point(i); q2_via_point(i)]));  
        first_itteration_for_loop = 0;
    else
        delta_x_hat(:,:,t_ind) = delta_x_hat(:,:,t_ind-1) + delta_T*((A(:,:,i)-F(:,:,i)*C(:,:,i))*delta_x_hat(:,:,t_ind-1) + B(:,:,i)*(u_input(:,:,t_ind-1) - u_op_each_itteration(:,:,t_ind-1)) + F(:,:,i)*(q - [q1_via_point(i); q2_via_point(i)]));
        delta_x_hat(1,:,t_ind) = q(1)-x_op(1,:,i);
        delta_x_hat(3,:,t_ind) = q(2)-x_op(3,:,i);
    end
    
    u_input(:,:,t_ind) = u_op(:,:,i) - K(:,:,i)*delta_x_hat(:,:,t_ind);
    [tout,qout] = ode45(@(time,x)simulatorofficial(time,x,u_input(:,:,t_ind),l1,l2,m1,m2,g,c1,c2),[t t+0.001],qout(end,:));
    q=qout(end,[1,3])';
    q1_visualize(t_ind+1) = q(1); % need t_ind+1 since "first_itteration_no_loop" takes first index position
    q2_visualize(t_ind+1) = q(2);
    q1d_visualize(t_ind+1) = qout(end,2);
    q2d_visualize(t_ind+1) = qout(end,4);
    q1_via_visualize(t_ind+1) = q1_via_point(i);
    q2_via_visualize(t_ind+1) = q2_via_point(i);
    
    x_current = l1*cos(q(1))+l2*cos(q(1)+q(2));
    y_current = l1*sin(q(1))+l2*sin(q(1)+q(1));
    position_current = [x_current; y_current];

    x_desired = l1*cos(q1_via_point(i))+l2*cos(q1_via_point(i)+q2_via_point(i));
    y_desired = l1*sin(q1_via_point(i))+l2*sin(q1_via_point(i)+q2_via_point(i));
    position_desired = [x_desired; y_desired];

    dist = norm(position_current - position_desired);

    % case 1: our desired point is a critical via point (i.e. point A, B, C, or D) and we are within 4mm of it, so we incriment holding_time
    if ((dist < required_dist) && ((i == lin_space_points - 1) || (i == lin_space_points*2 - 2) || (i == lin_space_points*3 - 3) || (i == lin_space_points*4 - 4)))
        holding_time = holding_time + 0.001;
    % case 2: our desired point is a critical via point (i.e. point A, B, C, or D) and we are not within 4mm of it, so we reset holding_time  
    elseif ((dist >= required_dist) && ((i == lin_space_points - 1) || (i == lin_space_points*2 - 2) || (i == lin_space_points*3 - 3) || (i == lin_space_points*4 - 4)))
        holding_time = 0;
    % case 3: our desired point is not a critical via point and we are within 4mm of it, so we can move to next point
    elseif dist < required_dist
        holding_time = 0;
        i = i + 1;
        ready_for_next_via_point = 1;
    end

    % we have held our position within 4mm for 5ms at a critical via point, so we can move to next point
    if holding_time == 0.005
        holding_time = 0;
        i = i + 1;
        ready_for_next_via_point = 1;
    end
    
    if i == lin_space_points*4 - 3
        break;
    end
    
    sprintf('t value = %d, i value = %d, dist = %d\n', t, i, dist)
    
end

params = [m1 m2 l1 l2 c1 c2];
visualize( params, tspan.', real(q1_visualize), real(q2_visualize), 'vis_plot.gif');

% figure; 
% subplot(1,2,1);
% plot(tspan, [real(q1_visualize), real(q1d_visualize), real(q2_visualize), real(q2d_visualize)]);
% legend('q1', 'q1d', 'q2', 'q2d');
% subplot(1,2,2);
% plot(tspan, [real(q1_via_visualize), real(q2_via_visualize)]);
% legend('q1via', 'q2via');