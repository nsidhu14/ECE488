clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

A = [0 0 -1;
    -1 0 0;
    3 -1 -3];

B = [0;
    0;
    1];

C = [0 1 0];

ev_desired = [-2 -2.0001 -2.00001];

F_transpose = place(A.', C.', ev_desired);

F = F_transpose.';

A_ss = [(A-F*C) F*C;
        zeros(3) A];

B_ss = [B;
        B];

C_ss = [zeros(1,3) C];

D_ss = [0];                           % Note: D_ss must have the same number of rows as C_ss

% Three sets of initial conditions x0 = [xhat1 xhat2 xhat3 x1 x2 x3]
x0_1 = [0; 0; 0; 1; 10; 100];
x0_2 = [0; 0; 0; 0.1; 0.01; 0.001];
x0_3 = [0; 0; 0; 5; 50; 50];

sys = ss(A_ss, B_ss, C_ss, D_ss);

tspan = 0:0.01:20;
u = zeros(size(tspan,2), 1);          % Note: u is a row vector of zeros (given in question) for every instance of t

% Three sets of simulations based on the initial conditons above
[y1, t1, x1] = lsim(sys, u, tspan, x0_1);
[y2, t2, x2] = lsim(sys, u, tspan, x0_2);
[y3, t3, x3] = lsim(sys, u, tspan, x0_3);

figure
plot(t1,x1);
title('Plot with first set of x0');
legend('xhat1', 'xhat2','xhat3','x1','x2','x3');

figure
plot(t2,x2);
title('Plot with second set of x0');
legend('xhat1', 'xhat2','xhat3','x1','x2','x3');

figure
plot(t3,x3);
title('Plot with third set of x0');
legend('xhat1', 'xhat2','xhat3','x1','x2','x3');


    