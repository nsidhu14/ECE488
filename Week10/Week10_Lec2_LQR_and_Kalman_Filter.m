clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

%% Deliverable 1 - using lqr to find K gain

A = [0 1 0 0;
    0 0 1 0;
    -3 1 2 3;
    2 1 0 0];

B = [0 0;
    0 0;
    1 2;
    0 2];

C = [1 0 0 0;
    0 0 0 1];

D = [0 0;                                       
    0 0];

% equal weighting on all states and inputs
Q1 = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

R1 = [1 0;
    0 1];

% minimizing u1 by increasing the weighting on R11
Q2 = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

R2 = [10 0;
    0 1];

% minimizing x4 by increasing the weighting on Q44
Q3 = [1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 100];

R3 = [1 0;
    0 1];

[K1, P1, ev1] = lqr(A,B,Q1,R1);
[K2, P2, ev2] = lqr(A,B,Q2,R2);
[K3, P3, ev3] = lqr(A,B,Q3,R3);

A_sf1 = A-B*K1;
B_sf1 = B;
C_sf1 = C-D*K1;
D_sf1 = D;

A_sf2 = A-B*K2; 
B_sf2 = B;
C_sf2 = C-D*K2;
D_sf2 = D;

A_sf3 = A-B*K3; 
B_sf3 = B;
C_sf3 = C-D*K3;
D_sf3 = D;

sys_sf1 = ss(A_sf1, B_sf1, C_sf1, D_sf1);
sys_sf2 = ss(A_sf2, B_sf2, C_sf2, D_sf2);
sys_sf3 = ss(A_sf3, B_sf3, C_sf3, D_sf3);

x0_sf = [1 1 1 1];

tspan = 0:0.01:10;
u = [zeros(size(tspan,2),1) zeros(size(tspan,2),1)];    % zero input
                                                                                                                                                        
[y_sf1, t_sf1, x_sf1] = lsim(sys_sf1, u, tspan, x0_sf);
[y_sf2, t_sf2, x_sf2] = lsim(sys_sf2, u, tspan, x0_sf);
[y_sf3, t_sf3, x_sf3] = lsim(sys_sf3, u, tspan, x0_sf);

figure
subplot(2,1,1);
plot(t_sf1, x_sf1);
title ('Plot of states using Q1 and R1');
legend('x1', 'x2', 'x3', 'x4');
subplot(2,1,2);
plot(t_sf1, -K1*x_sf1');
title ('Plot of inputs using Q1 and R1');
legend('u1', 'u2');

figure
subplot(2,1,1);
plot(t_sf2, x_sf2);
title ('Plot of states using Q2 and R2');
legend('x1', 'x2', 'x3', 'x4');
subplot(2,1,2);
plot(t_sf2, -K2*x_sf2');
title ('Plot of inputs using Q2 and R2');
legend('u1', 'u2');

figure
subplot(2,1,1);
plot(t_sf3, x_sf3);
title ('Plot of states using Q3 and R3');
legend('x1', 'x2', 'x3', 'x4');
subplot(2,1,2);
plot(t_sf3, -K3*x_sf3');
title ('Plot of inputs using Q3 and R3');
legend('u1', 'u2');

%% Deliverable 2 - using Kalman filter to find F gain

a = 0.1; % std dev of both 'n' and 'w'
b = 0; % mean of both 'n' and 'w'
n = a.*randn(4,1001) + b; % n has 4 rows since its part of 'x_dot' which has 4 states
w = a.*randn(2,1001) + b; % w has 2 rows since its part of 'y' which has 2 outputs

Q_kf = eye(4)*a^2; % diagonal matrix containing variance of 'n' as the diagonal entries
R_kf = eye(2)*a^2; % diagonal matrix containing variance of 'w' as the diagonal entries

K_kf = K1; % using the gain 'K' obtained from the first case in deliverable 1

[F_transpose_kf, P_kf, ev_kf] = lqr(A.',C.',Q_kf,R_kf); % Note: this returns F transpose
F_kf = transpose(F_transpose_kf);

A_kf = [A -B*K_kf;
    F_kf*C A-F_kf*C-B*K_kf];

B_kf = [B eye(4) zeros(4,2);
    B zeros(4) F_kf];

C_kf = [C zeros(2,4)];

D_kf = [zeros(2) zeros(2,4) eye(2)];

sys_kf = ss(A_kf, B_kf, C_kf, D_kf);

x0_kf = [1 1 1 1 1 1 1 1];

tspan_kf = 0:0.01:10;
u_kf = [zeros(size(tspan_kf,2),1) zeros(size(tspan_kf,2),1) n.' w.'];   
                                                                                                                                                        
[y_kf, t_kf, x_kf] = lsim(sys_kf, u_kf, tspan_kf, x0_kf);

figure
subplot(2,1,1);
plot(t_kf, x_kf);
title ('Plot of states using Q\_kf and R\_kf');
legend('x1', 'x2', 'x3', 'x4', 'x1hat', 'x2hat', 'x3hat', 'x4hat');
subplot(2,1,2);
plot(t_kf, -K_kf*x_kf(:,(1:4))');
title ('Plot of inputs using Q\_kf and R\_kf');
legend('u1', 'u2');

