clear all; 
close all; 

%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

%% State Space Linearization

syms theta theta_dot theta_d_dot tau 
eqn = theta_d_dot + theta_dot + sin(theta) == tau;
theta_d_dot_val = solve(eqn, theta_d_dot);

x = [theta; theta_dot];
x_dot = [theta_dot; theta_d_dot_val];
y = theta;

A = jacobian(x_dot,x);
B = jacobian(x_dot,tau);
C = jacobian(y,x);
D = jacobian(y,tau);

% Substituing the equilibrium points and converting from syms to numbers type 
A = double(subs(A, [theta,theta_dot,tau], [0,0,0]));
B = double(subs(B, [theta,theta_dot,tau], [0,0,0]));
C = double(subs(C, [theta,theta_dot,tau], [0,0,0]));
D = double(subs(D, [theta,theta_dot,tau], [0,0,0]));

tspan = [0:0.01:10];
x0 = [0;0];

% Plotting the linearized state space output
[t,x_output_linear] = ode45(@(t,x) Week2_Lec2_Func_ss_sys_linear(t,x,A,B,tau),tspan,x0); 
figure
plot(t,x_output_linear(:,1))
title('Linearized Output');

% Plotting the non-linear state space output
[t,x_output_nonlinear] = ode45(@(t,x) Week2_Lec2_Func_ss_sys_nonlinear(t,x,theta,theta_dot,tau,x_dot),tspan,x0); 
figure
plot(t,x_output_nonlinear(:,1))
title('Non-Linearized Output');