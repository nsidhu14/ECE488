clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

s = tf('s');

% start with controllable cononical form
A_bar = [0 1 0;
         0 0 1;
        -1 -3 -3];

B_bar = [0;
         0;
         1];
     
C_bar = [1 0 0];

D_bar = [0];

K_bar = [7 9 3];

ev_original = eig(A_bar);          % the original eigenvalues are all at -1.

ev_final = eig(A_bar-B_bar*K_bar); % the final eigenvalues are all at -2.

P = [0 1 0;                        % arbitrary transformation matrix
    1 1 0;
    0 0 1];

P_inv = inv(P);

A = P_inv*A_bar*P;
B = P_inv*B_bar;
K = K_bar*P;
C = C_bar*P;
D = D_bar;

