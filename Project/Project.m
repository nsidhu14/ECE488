clear all; 
close all;

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
c1 = 1;
c2 = 1;
g = 1;

syms q1 q2 q1_d q2_d q1_dd q2_dd tau1 tau2
eqns = [tau1 == ((m1*l1^2)/3 + (m2*l2^2)/12 + m2*(l1^2+(l2^2)/4+l1*l2*cos(q2)))*q1_dd + ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q2_dd - m2*l1*l2*sin(q2)*q1_d*q2_d - (m2*l1*l2*sin(q2))/2*q2_d^2 + ((m1*l1)/2+m2*l1)*g*cos(q1) + (m2*l2)/2*g*cos(q1+q2) + c1*q1_d , tau2 == ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q1_dd + (m2*l2^2)/3*q2_dd + (m2*l1*l2*sin(q2))/2*q1_d^2 + (m2*l2)/2*g*cos(q1+q2) + c2*q2_d];
dd_vals = solve(eqns, [q1_dd q2_dd]);

x = [q1; q1_d; q2; q2_d];
x_dot = [q1_d; dd_vals.q1_dd; q2_d; dd_vals.q2_dd];
u = [tau1; tau2];
y = [q1; q2];

A = jacobian(x_dot,x);
B = jacobian(x_dot,u);
C = jacobian(y,x);
D = jacobian(y,u);

x_eq = [pi/2; 0; 0; 0];
u_eq = [0; 0];

x0 = [pi/2;0;0;0];

% Substituing the equilibrium points and converting from syms to numbers type 
A = double(subs(A, [x.' u.'], [x_eq.' u_eq.']));       % Note: use .' to transpose x and u into one long row vector
B = double(subs(B, [x.' u.'], [x_eq.' u_eq.']));
C = double(subs(C, [x.' u.'], [x_eq.' u_eq.']));
D = double(subs(D, [x.' u.'], [x_eq.' u_eq.']));

delta_T = 0.5;
time_limit = 20;
tspan = [0:delta_T:time_limit];

% Plotting the non-linear state space output
[t,x_output_nonlinear] = ode45(@(t,x) Func_ss_sys_nonlinear_no_syms(t,x),tspan,x0);
figure
plot(t,x_output_nonlinear)
title('Non-Linearized Output');
legend('q1','q1dot','q2','q2dot');
xlim([0 time_limit]);


% Plotting the linearized state space output
% Note: The linearized output from ode45 gives the delta_x vector, yet x = x_eq + delta_x
[t,delta_x_output_linear] = ode45(@(t,x) Func_ss_sys_linear(t,x,A,B,u),tspan,x0);
x_output_linear = [x_eq(1)+delta_x_output_linear(:,1), x_eq(2)+delta_x_output_linear(:,2), x_eq(3)+delta_x_output_linear(:,3), x_eq(4)+delta_x_output_linear(:,4)];
figure
plot(t,x_output_linear)
title('Linearized Output');
legend('q1','q1dot','q2','q2dot');
xlim([0 time_limit]);

% Pole/Zero Map for linearized state-space system
ss_sys = ss(A,B,C,D);     % Creating the state-space model
figure
pzmap(ss_sys);
title('Pole-Zero Map of Linearized State-Space System');

%% Creating PD controller (underdamped)

% Plant transfer function (this is a 2x2 t.f.)
P = tf(ss_sys); 

figure
pzmap(P);
title('Pole-Zero Map of Linearized Plant');


%% State Feedback
ev_desired = [-1 -1 -1 -1];
ev_A = eig(A);

for i = 1:size(ev_A,1)
    if ev_A(i) < 0
        ev_desired(i) = ev_A(i);
    end
end

K = place(A, B, ev_desired);

U0 = [0 0];
u_input = U0 -K*(x-x0);

%%
% s = tf('s');
% 
% Kd11 = 1;
% Kp11 = 1;
% 
% Kd12 = 1;
% Kp12 = 1;
% 
% Kd21 = 1;
% Kp21 = 1;
% 
% Kd22 = 1;
% Kp22 = 1;
% 
% C11 = Kd11*s + Kp11;
% C12 = Kd12*s + Kp12;
% C21 = Kd21*s + Kp21;
% C22 = Kd22*s + Kp22;
% 
% C = [C11 C12;
%      C21 C22];
%  
% closed_loop_PD_tf = feedback(P*C,eye(2));
% 
% closed_loop_PD_ss_sys = ss(closed_loop_PD_tf);
% 
% [A_PD, B_PD, C_PD, D_PD] = ssdata(closed_loop_PD_ss_sys);
% 
% % Plotting the linearized PD state space output
% % Note: The linearized output from ode45 gives the delta_x vector, yet x = x_eq + delta_x
% [t_PD,delta_x_output_linear_PD] = ode45(@(t,x) Func_ss_sys_linear(t,x,A_PD,B_PD,u),tspan,x0);
% x_output_linear_PD = [x_eq(1)+delta_x_output_linear_PD(:,1), x_eq(2)+delta_x_output_linear_PD(:,2), x_eq(3)+delta_x_output_linear_PD(:,3), x_eq(4)+delta_x_output_linear_PD(:,4)];
% figure
% plot(t_PD,x_output_linear_PD)
% title('Linearized PD Output');
% legend('q1','q1dot','q2','q2dot');
% xlim([0 time_limit]);
