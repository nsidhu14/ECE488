%% Function to create the non-linear state-space system

function dxdt = Week3_Tut_Func_ss_sys_nonlinear_no_syms(t,x)

   u1 = 0;
   u2 = 0;
   
   q1 = x(1);
   q1_d = x(2);
   q2 = x(3);
   q2_d = x(4);
   tau1 = u1;
   tau2 = u2;
   
   x1_dot = double(q1_d);
   x2_dot = double(-(3*(4*q2_d - 4*q1_d + 4*tau1 - 4*tau2 - 6*cos(q1) + 2*q1_d^2*sin(q2) + 2*q2_d^2*sin(q2) + 6*q2_d*cos(q2) - 6*tau2*cos(q2) + 3*cos(q1 + q2)*cos(q2) + 4*q1_d*q2_d*sin(q2) + 3*q1_d^2*cos(q2)*sin(q2)))/(9*cos(q2)^2 - 16));
   x3_dot = double(q2_d);
   x4_dot = double((3*(20*q2_d - 4*q1_d + 4*tau1 - 20*tau2 + 8*cos(q1 + q2) - 6*cos(q1) - 9*cos(q1)*cos(q2) + 10*q1_d^2*sin(q2) + 2*q2_d^2*sin(q2) - 6*q1_d*cos(q2) + 12*q2_d*cos(q2) + 6*tau1*cos(q2) - 12*tau2*cos(q2) + 3*cos(q1 + q2)*cos(q2) + 4*q1_d*q2_d*sin(q2) + 6*q1_d^2*cos(q2)*sin(q2) + 3*q2_d^2*cos(q2)*sin(q2) + 6*q1_d*q2_d*cos(q2)*sin(q2)))/(9*cos(q2)^2 - 16));

%    x1_dot = double(subs(x_dot(1), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
%    x2_dot = double(subs(x_dot(2), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
%    x3_dot = double(subs(x_dot(3), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));
%    x4_dot = double(subs(x_dot(4), [q1,q1_d,q2,q2_d,tau1,tau2], [x(1),x(2),x(3),x(4),u1,u2]));

   dxdt = [x1_dot; x2_dot; x3_dot; x4_dot];
end