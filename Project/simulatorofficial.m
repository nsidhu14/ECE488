%% Function to create the non-linear state-space system without using symbolic substitution

function dxdt = simulatorofficial(t,x,U,l1,l2,m1,m2,g,c1,c2)
   u1 = U(1); % Set the input1 value
   u2 = U(2); % Set the input2 value
   tau1 = u1;
   tau2 = u2;
   
   q1 = x(1);
   q1_d = x(2);
   q2 = x(3);
   q2_d = x(4);
   
   x1_dot = double(q1_d);
   x2_dot = double((3*(4*l2*tau1 - 4*l2*tau2 - 6*l1*tau2*cos(q2) - 4*c1*l2*q1_d + 4*c2*l2*q2_d + 6*c2*l1*q2_d*cos(q2) + 2*l1*l2^2*m2*q1_d^2*sin(q2) + 2*l1*l2^2*m2*q2_d^2*sin(q2) - 2*g*l1*l2*m1*cos(q1) - 4*g*l1*l2*m2*cos(q1) + 3*g*l1*l2*m2*cos(q1 + q2)*cos(q2) + 3*l1^2*l2*m2*q1_d^2*cos(q2)*sin(q2) + 4*l1*l2^2*m2*q1_d*q2_d*sin(q2)))/(4*l1^2*l2*m1 + 12*l1^2*l2*m2 - 9*l1^2*l2*m2*cos(q2)^2));
   x3_dot = double(q2_d);
   x4_dot = double(-(3*(4*l2^2*m2*tau1 - 12*l1^2*m2*tau2 - 4*l1^2*m1*tau2 - 4*l2^2*m2*tau2 - 4*c1*l2^2*m2*q1_d + 4*c2*l1^2*m1*q2_d + 12*c2*l1^2*m2*q2_d + 4*c2*l2^2*m2*q2_d + 2*l1*l2^3*m2^2*q1_d^2*sin(q2) + 6*l1^3*l2*m2^2*q1_d^2*sin(q2) + 2*l1*l2^3*m2^2*q2_d^2*sin(q2) + 6*g*l1^2*l2*m2^2*cos(q1 + q2) - 4*g*l1*l2^2*m2^2*cos(q1) + 6*l1*l2*m2*tau1*cos(q2) - 12*l1*l2*m2*tau2*cos(q2) + 2*l1^3*l2*m1*m2*q1_d^2*sin(q2) + 4*l1*l2^3*m2^2*q1_d*q2_d*sin(q2) + 3*g*l1*l2^2*m2^2*cos(q1 + q2)*cos(q2) + 6*l1^2*l2^2*m2^2*q1_d^2*cos(q2)*sin(q2) + 3*l1^2*l2^2*m2^2*q2_d^2*cos(q2)*sin(q2) - 6*c1*l1*l2*m2*q1_d*cos(q2) + 12*c2*l1*l2*m2*q2_d*cos(q2) - 6*g*l1^2*l2*m2^2*cos(q1)*cos(q2) + 2*g*l1^2*l2*m1*m2*cos(q1 + q2) - 2*g*l1*l2^2*m1*m2*cos(q1) + 6*l1^2*l2^2*m2^2*q1_d*q2_d*cos(q2)*sin(q2) - 3*g*l1^2*l2*m1*m2*cos(q1)*cos(q2)))/(12*l1^2*l2^2*m2^2 + 4*l1^2*l2^2*m1*m2 - 9*l1^2*l2^2*m2^2*cos(q2)^2));

   dxdt = [x1_dot; x2_dot; x3_dot; x4_dot];
end