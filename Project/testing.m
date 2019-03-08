clear all; 
close all;

m1 = 0.375;
m2 = 0.375;
l1 = 0.15;
l2 = 0.15;
c1 = 2;
c2 = 2;
g = 3.7;


x = [0.1 0.2 0.2 0.1];
y = [0.2 0.2 0.1 0.1];

q1=zeros(1,4);
q2=zeros(1,4);
for i=1:4
    q1(i) = atan2(y(i),x(i))-acos((x(i)^2 + y(i)^2 + l1^2 - l2^2)/(2*l1*sqrt(x(i)^2+y(i)^2)));
    q2(i) = acos((x(i)^2+y(i)^2-l1^2-l2^2)/(2*l1*l2));
end

%% Total angular displacement

TAD = abs(q1(4)-q1(3))+abs(q1(3)-q1(2))+abs(q1(2)-q1(1))+abs(q1(1)-q1(4))+ abs(q2(4)-q2(3))+abs(q2(3)-q2(2))+abs(q2(2)-q2(1))+abs(q2(1)-q2(4))

%% EOAT Position plot
q1_plot(1,:) = linspace(q1(1),q1(2));
q1_plot(2,:) = linspace(q1(2),q1(3));
q1_plot(3,:) = linspace(q1(3),q1(4));
q1_plot(4,:) = linspace(q1(4),q1(1));

q2_plot(1,:) = linspace(q2(1),q2(2));
q2_plot(2,:) = linspace(q2(2),q2(3));
q2_plot(3,:) = linspace(q2(3),q2(4));
q2_plot(4,:) = linspace(q2(4),q2(1));

x_plot = [l1*cos(q1_plot(1,:))+l2*cos(q1_plot(1,:)+q2_plot(1,:)),l1*cos(q1_plot(2,:))+l2*cos(q1_plot(2,:)+q2_plot(2,:)),l1*cos(q1_plot(3,:))+l2*cos(q1_plot(3,:)+q2_plot(3,:)),l1*cos(q1_plot(4,:))+l2*cos(q1_plot(4,:)+q2_plot(4,:))];
y_plot = [l1*sin(q1_plot(1,:))+l2*sin(q1_plot(1,:)+q2_plot(1,:)),l1*sin(q1_plot(2,:))+l2*sin(q1_plot(2,:)+q2_plot(2,:)),l1*sin(q1_plot(3,:))+l2*sin(q1_plot(3,:)+q2_plot(3,:)),l1*sin(q1_plot(4,:))+l2*sin(q1_plot(4,:)+q2_plot(4,:))];

figure
hold on
scatter(x_plot,y_plot);
plot(0,0,'r*');

% a = l1;
% b = l2;
% c = sqrt(5)/10;
% 
% C = acos((a^2 + b^2 - c^2)/(2*a*b));
% A = acos((b^2 + c^2 - a^2)/(2*b*c));
% B = acos((c^2 + a^2 - b^2)/(2*c*a));
% 
% B_max = asin(0.20/c);
% 
% q1_start = B_max - B;
% q2_start = pi - C;
% 
% x0 = 0;
% y0 = 0;
% x1 = a*cos(q1_start);
% y1 = a*sin(q1_start);
% x2 = c*cos(B_max);
% y2 = c*sin(B_max);
% 
% figure 
% hold on 
% plot(x0,y0,'r*');
% plot(x1,y1,'r*');
% plot(x2,y2,'r*');
% plot([x0 x1], [y0 y1]);
% plot([x1 x2], [y1 y2]);
% hold off
% 
% 

q1_start = q1_plot(1,1);
q2_start = q1_plot(1,1);

syms q1 q2 q1_d q2_d q1_dd q2_dd tau1 tau2
eqns = [tau1 == ((m1*l1^2)/3 + (m2*l2^2)/12 + m2*(l1^2+(l2^2)/4+l1*l2*cos(q2)))*q1_dd + ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q2_dd - m2*l1*l2*sin(q2)*q1_d*q2_d - (m2*l1*l2*sin(q2))/2*q2_d^2 + ((m1*l1)/2+m2*l1)*g*cos(q1) + (m2*l2)/2*g*cos(q1+q2) + c1*q1_d , tau2 == ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q1_dd + (m2*l2^2)/3*q2_dd + (m2*l1*l2*sin(q2))/2*q1_d^2 + (m2*l2)/2*g*cos(q1+q2) + c2*q2_d];
dd_vals = solve(eqns, [q1_dd q2_dd]);
tau_vals = solve(eqns, [tau1 tau2]);

x = [q1; q1_d; q2; q2_d];
x_dot = [q1_d; dd_vals.q1_dd; q2_d; dd_vals.q2_dd];
u = [tau1; tau2];
y = [q1; q2];

A = jacobian(x_dot,x);
B = jacobian(x_dot,u);
C = jacobian(y,x);
D = jacobian(y,u);

q1_next = q1_plot(1,2);
q2_next = q2_plot(1,2);

x_op = [q1_next; 0; q2_next; 0];
tau1_op = double(subs(tau_vals.tau1, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_next 0 0 q2_next 0 0]));
tau2_op = double(subs(tau_vals.tau2, [q1 q1_d q1_dd q2 q2_d q2_dd ], [q1_next 0 0 q2_next 0 0]));
u_op = [tau1_op; tau2_op];

x0 = [q1_start;0;q2_start;0];

% Substituing the equilibrium points and converting from syms to numbers type 
A = double(subs(A, [x.' u.'], [x_op.' u_op.']));       % Note: use .' to transpose x and u into one long row vector
B = double(subs(B, [x.' u.'], [x_op.' u_op.']));
C = double(subs(C, [x.' u.'], [x_op.' u_op.']));
D = double(subs(D, [x.' u.'], [x_op.' u_op.']));

% %% State Feedback
ev_desired = [-1 -1 -1 -1];
ev_A = eig(A);

for i = 1:size(ev_A,1)
    if ev_A(i) < 0
        ev_desired(i) = ev_A(i);
    end
end

K = place(A, B, ev_desired);

u_input = u_op -K*(x-x_op);

delta_T = 0.5;
time_limit = 20;
tspan = [0:delta_T:time_limit];

% Plotting the non-linear state space output
[t,x_output_nonlinear] = ode45(@(t,x)simulatorofficial(t,x,u_input,l1,l2,m1,m2,g,c1,c2),tspan,x0);
figure
plot(t,x_output_nonlinear)
title('Non-Linearized Output');
legend('q1','q1dot','q2','q2dot');
xlim([0 time_limit]);

% 
% % Plotting the linearized state space output
% % Note: The linearized output from ode45 gives the delta_x vector, yet x = x_eq + delta_x
% [t,delta_x_output_linear] = ode45(@(t,x) Func_ss_sys_linear(t,x,A,B,u),tspan,x0);
% x_output_linear = [x_eq(1)+delta_x_output_linear(:,1), x_eq(2)+delta_x_output_linear(:,2), x_eq(3)+delta_x_output_linear(:,3), x_eq(4)+delta_x_output_linear(:,4)];
% figure
% plot(t,x_output_linear)
% title('Linearized Output');
% legend('q1','q1dot','q2','q2dot');
% xlim([0 time_limit]);
% 
% % Pole/Zero Map for linearized state-space system
% ss_sys = ss(A,B,C,D);     % Creating the state-space model
% figure
% pzmap(ss_sys);
% title('Pole-Zero Map of Linearized State-Space System');
% 
% %% Creating PD controller (underdamped)
% 
% % Plant transfer function (this is a 2x2 t.f.)
% P = tf(ss_sys); 
% 
% figure
% pzmap(P);
% title('Pole-Zero Map of Linearized Plant');
% 
% 


