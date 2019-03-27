clear all; 
close all;

% To do:
% - will need to have non-zero velocity at each via point
% - stuck here : simulate y_output for each itteration
% - replace manual pole placement with LQR / Kalman
% - why is our system unstable?
% - why does our simulation take so long?
% - how to use visualize?
% - time limit?

m1 = 0.375;
m2 = 0.375;
l1 = 0.15;
l2 = 0.15;
c1 = 2;
c2 = 2;
g = 3.7;

x_boundary = [0 0 0.22 0.22; 0 0.22 0.22 0]; 
y_boundary = [0 0.22 0.22 0; 0.22 0.22 0 0]; 

lin_space_points = 5;

delta_T = 0.001;
t_max = 10;
tspan = 0:delta_T:t_max;

% initializing all variables for each via point
A = sym(zeros(4,4,lin_space_points*4-3));
B = sym(zeros(4,2,lin_space_points*4-3));
C = sym(zeros(2,4,lin_space_points*4-3));
D = sym(zeros(2,2,lin_space_points*4-3));

ev_A = zeros(4,1,lin_space_points*4-3);
ev_A = complex(ev_A);

K = zeros(2,4,lin_space_points*4-3);
K = complex(K);

F = zeros(4,2,lin_space_points*4-3);
F = complex(F);

x_op = zeros(4,1,lin_space_points*4-3);
tau1_op = zeros(lin_space_points*4-3,1);
tau2_op = zeros(lin_space_points*4-3,1);
u_op = zeros(2,1,lin_space_points*4-3);
x0 = zeros(4,1,lin_space_points*4-3);

y_output = zeros(2,1,length(tspan));
x_hat = zeros(4,1,length(tspan));
u_input = zeros(2,1,length(tspan));

q1_visualize = zeros(size(tspan,2),1);
q2_visualize = zeros(size(tspan,2),1);

Q_lqr = zeros(4,4,lin_space_points*4-3);
R_lqr = zeros(2,2,lin_space_points*4-3);
Q_kf = zeros(4,4,lin_space_points*4-3);
R_kf = zeros(2,2,lin_space_points*4-3);

mean = 0;
std_dev = 1/3;
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
q1_via_point = [q1_via_point(1,(1:end-1)) q1_via_point(2,(1:end-1)) q1_via_point(3,(1:end-1)) q1_via_point(4,:)];
q2_via_point = [q2_via_point(1,(1:end-1)) q2_via_point(2,(1:end-1)) q2_via_point(3,(1:end-1)) q2_via_point(4,:)];

%% Plotting

figure
hold on
scatter(x_plot_via_point,y_plot_via_point, '.b');
plot(x_boundary, y_boundary, '-r');
plot(x_coord, y_coord, '+g');

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
i_current = i; i_prev = 0;

for t=0.001:0.001:t_max
    t_ind = round(t/delta_T);
    if i_current ~= i_prev
        A(:,:,i) = jacobian(x_dot,x);
        B(:,:,i) = jacobian(x_dot,u);
        C(:,:,i) = jacobian(y,x);
        D(:,:,i) = jacobian(y,u);

        x_op(:,:,i) = [q1_via_point(i+1); 0; q2_via_point(i+1); 0];
        tau1_op(i) = double(subs(tau_vals.tau1, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(i+1) 0 0 q2_via_point(i+1) 0 0]));
        tau2_op(i) = double(subs(tau_vals.tau2, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(i+1) 0 0 q2_via_point(i+1) 0 0]));
        u_op(:,:,i) = [tau1_op(i); tau2_op(i)];

        x0(:,:,i) = [q1_via_point(i); 0; q2_via_point(i); 0];

        % Substituing the equilibrium points and converting from syms to numbers type 
        A(:,:,i) = double(subs(A(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));       
        B(:,:,i) = double(subs(B(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
        C(:,:,i) = double(subs(C(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
        D(:,:,i) = double(subs(D(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));

        x_next = l1*cos(q1_via_point(i+1))+l2*cos(q1_via_point(i+1)+q2_via_point(i+1));
        y_next = l1*sin(q1_via_point(i+1))+l2*sin(q1_via_point(i+1)+q2_via_point(i+1));
        plot(x_next, y_next, 'om');

%         ev_A(:,:,i) = eig(A(:,:,i));
%         ev_desired_K = [-10 -20 -30 -40];
%         ev_desired_F = [-20 -40 -60 -80];
% 
%         for j = 1:size(ev_A(:,:,i),1)
%             if ev_A(j,:,i) < 0
%                 ev_desired_K(j) = ev_A(j,:,i);
%             end
%             ev_desired_F(j) = 2*ev_desired_K(j);
%         end
% 
%         K(:,:,i) = place(double(A(:,:,i)), double(B(:,:,i)), ev_desired_K);
%         F(:,:,i) = transpose(place(double(transpose(A(:,:,i))), double(transpose(C(:,:,i))), ev_desired_F));
       
        Q_lqr(:,:,i) = [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 1];
        R_lqr(:,:,i) = [1 0;
            0 1];
        [K_lqr, P_lqr, ev_lqr] = lqr(double(A(:,:,i)), double(B(:,:,i)), Q_lqr(:,:,i), R_lqr(:,:,i));
        K(:,:,i) = K_lqr;
        
        Q_kf(:,:,i) = eye(4)*variance;
        R_kf(:,:,i) = eye(2)*variance;
        [F_transpose_kf, P_kf, ev_kf] = lqr(double(A(:,:,i)).', double(C(:,:,i)).', Q_kf(:,:,i), R_kf(:,:,i)); % Note: this returns F transpose
        F(:,:,i) = transpose(F_transpose_kf);
    end
    i_prev = i_current;
    % accounting for first itteration
    if i == 1
        x_hat(:,:,t_ind) = x0(:,:,i);
        y_output(:,:,t_ind) = [q1_via_point(i); q2_via_point(i)];
        x_hat(:,:,t_ind+1) = x_hat(:,:,t_ind) + delta_T*((A(:,:,i)-F(:,:,i)*C(:,:,i)-B(:,:,i)*K(:,:,i))*x_hat(:,:,t_ind) + F(:,:,i)*y_output(:,:,t_ind));
        u_input(:,:,t_ind+1) = u_op(:,:,i) - K(:,:,i)*(x_hat(:,:,t_ind+1)-x_op(:,:,i));
        q1_visualize(t_ind) = q1_via_point(i);
        q2_visualize(t_ind) = q2_via_point(i);
        i = i + 1;
        i_current = i_current + 1;
    else
        % finding y_output(i) using output of non-linear ode45, 
        % subbing u_input(:,:,i) which is already found since we always solve u_input for next itteration 
        [t_nl,x_nl] = ode45(@(time,x)simulatorofficial(time,x,u_input(:,:,t_ind),l1,l2,m1,m2,g,c1,c2),[t t+0.001],x0(:,:,i)); 
        y_output(:,:,t_ind) = x_nl(end,[1,3])';
        q1_visualize(t_ind) = x_nl(end,1)';
        q2_visualize(t_ind) = x_nl(end,3)';
        x_hat(:,:,t_ind+1) = x_hat(:,:,t_ind) + delta_T*((A(:,:,i)-F(:,:,i)*C(:,:,i)-B(:,:,i)*K(:,:,i))*x_hat(:,:,t_ind) + F(:,:,i)*y_output(:,:,t_ind));
        u_input(:,:,t_ind+1) = u_op(:,:,i) - K(:,:,i)*(x_hat(:,:,t_ind+1)-x_op(:,:,i));
        
        if isequal((abs(y_output(:,:,t_ind)-[q1_via_point(i); q2_via_point(i)]) < [0.2; 0.2]),[1;1])
            i = i + 1;
            i_current = i_current + 1;
        end
    end
    if i == lin_space_points*4 - 3
        break;
    end
    sprintf('t value = %d, i value = %d, prev i = %d\n', t, i, i_prev)
end

params = [m1 m2 l1 l2 c1 c2];
visualize( params, tspan.', real(q1_visualize), real(q2_visualize), 'vis_plot.gif');