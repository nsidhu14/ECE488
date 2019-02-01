%% Function to create the non-linear state-space system

function dxdt = Week2_Lec2_Func_ss_sys_nonlinear(t,x,theta,theta_dot,tau,x_dot)
% pass in the symbolic terms 'theta', 'theta_dot', and 'tau' to replace
% them with the numbers of 'x(1)', 'x(2)', and 'u'

   u = 0.1;
   x1_dot = double(subs(x_dot(1), [theta,theta_dot,tau], [x(1),x(2),u]));
   x2_dot = double(subs(x_dot(2), [theta,theta_dot,tau], [x(1),x(2),u]));
   dxdt = [x1_dot; x2_dot];
end