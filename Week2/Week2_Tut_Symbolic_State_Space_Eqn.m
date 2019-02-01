clear all; 
close all; 

%% Solving the mass-spring-damper system symbolically and then using ode45

% defining the constants
M = 1;
B = 1;
K = 1;

syms d d_dot d_double_dot u 
eqn = M*d_double_dot + B*d_dot + K*d == u;
d_double_dot_val = solve(eqn, d_double_dot); % solved for d_double_dot with respect to "d" and "d_dot"

% solving for the 4 matrices, using double() to convert from syms type matrix to
% a double type matrix so it can be used in ode45() which does not accept
% symbolic inputs
x = [d; d_dot]
x_dot = [d_dot; d_double_dot_val];
y = d;

A = double(jacobian(x_dot,x));
B = double(jacobian(x_dot,u));
C = double(jacobian(y,x));
D = double(jacobian(y,u));

tspan = [0 20];
x0 = [0;0];

[t,x] = ode45(@(t,x) Week2_Tut_Func_ss_sys(t,x,A,B,u),tspan,x0); % first argument is the state space system x = A*x + B*u as a function

plot(t,x(:,1))
    