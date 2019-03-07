clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

%% Part 1: Combined system of state feedback and state estimator

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

D = [0 0;                                       % Note: D is 2x2 since we have 2 inputs and 2 outputs
    0 0];

ev_desired_K = [-1 -2 -3 -4];                   % state feedback ev's are chosen arbitrarily
ev_desired_F = [-2 -4 -6 -8];                   % state estimator ev's must be atleast 2x faster than state feedback ev's

K = place(A, B, ev_desired_K);
F_transpose = place(A.', C.', ev_desired_F);
F = transpose(F_transpose);

A_sf = A-B*K;                                   % 'A' matrix for state feedback only

A_combined = [A -B*K;
    F*C A-F*C-B*K];

B_combined = [B;
    B];
    
C_combined = [C -D*K];

D_combined = D;

sys_sf = ss(A_sf, B, C, D);
sys_combined = ss(A_combined, B_combined, C_combined, D_combined);

x0_sf = [0 0 0 0];                              % state feedback 'A' matirx is 4x4, thus need 4 IC's
x0_combined_all_zeros = [0 0 0 0 0 0 0 0];      % combined 'A' matrix is 8x8, thus need 8 IC's (4 for 'x' and 4 for 'x_hat')
x0_combined_diff_inputs = [1 2 3 4 5 6 7 8];

tspan = 0:0.01:20;
u = [ones(size(tspan,2),1) ones(size(tspan,2),1)];   % use ones() since 'u' is a step input & need 2 columns since 2 inputs 
                                                     % size(tspan,2) to get its number of columns which holds its length
                                                     % ones(size(tspan,2),1) use the ,1 so we only create a single column of ones instead of a square matrix
                                                                                                      
[y_sf, t_sf, x_sf] = lsim(sys_sf, u, tspan, x0_sf);
[y_combined_all_zeros, t_combined_all_zeros, x_combined_all_zeros] = lsim(sys_combined, u, tspan, x0_combined_all_zeros);
[y_combined_diff_inputs, t_combined_diff_inputs, x_combined_diff_inputs] = lsim(sys_combined, u, tspan, x0_combined_diff_inputs);

figure
plot(t_sf,x_sf);
title('Plot of state feedback without state estimator using all zero initial conditions');
legend('x1', 'x2', 'x3', 'x4');

figure
plot(t_combined_all_zeros,x_combined_all_zeros);
title('Plot of state feedback with state estimator using all zero initial conditions');
legend('x1', 'x2', 'x3', 'x4', 'xhat1', 'xhat2', 'xhat3', 'xhat4');

figure
plot(t_combined_diff_inputs,x_combined_diff_inputs);
title('Plot of state feedback with state estimator using different initial conditions');
legend('x1', 'x2', 'x3', 'x4', 'xhat1', 'xhat2', 'xhat3', 'xhat4');


%% Part 2: Regulator 

% z has states [x_dot;   which is a 6x1 matrix
%               n_dot]

A_aug = [A zeros(4,2);
    C zeros(2)];

B_aug = [B;
    zeros(2)];

C_aug = [C zeros(2)];

ev_desired_F_prime = [-1 -2 -3 -4 -5 -6];

F_prime = place(A_aug, B_aug, ev_desired_F_prime);

closed_loop_ev = eig(A_aug-B_aug*F_prime);  % verifying ev's

F1_prime = F_prime(:,(1:4));                % F1_prime is 2x4
F2_prime = F_prime(:,(5:6));                % F2_prime is 2x2 


% the final system has states [x;   which is a 6x1 matrix
%                              n];

A_system = [A-B*F1_prime -B*F2_prime;
    C zeros(2)];

B_system = [zeros(4,2);
    -eye(2)];

C_system = [C zeros(2)];

D_system = D;

sys_system = ss(A_system, B_system, C_system, D_system);

x0_system = [0 0 0 0 0 0];

tspan = 0:0.01:20;
u = [ones(size(tspan,2), 1) ones(size(tspan,2), 1)];

[y_system, t_system, x_system] = lsim(sys_system, u, tspan, x0_system);

figure
plot(t_system,x_system);
title('Plot of state feedback with integrator without estimator using all zero initial conditions');
legend('x1', 'x2', 'x3', 'x4', 'n1', 'n2');

