%% Function to create the non-linear state-space system

function dxdt = Week3_Tut_Func_ss_sys_nonlinear(t,x,q1,q1_d,q2,q2_d,tau1,tau2,x_dot)
% pass in the symbolic terms 'theta', 'theta_dot', and 'tau' to replace
% them with the numbers of 'x(1)', 'x(2)', and 'u'

   u1 = 1;
   u2 = 1;
   x1_dot = double(subs(x_dot(1), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
   x2_dot = double(subs(x_dot(2), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
   x3_dot = double(subs(x_dot(3), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
   x4_dot = double(subs(x_dot(4), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
   dxdt = [x1_dot; x2_dot; x3_dot; x4_dot];
end