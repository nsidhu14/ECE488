clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

%% Solving for P (similarity transformation matrix) given the A,B,C matrices

A = [0 0 -1;
    -1 -1 1;
    0 -1 -2];

B = [1;
    0;
    2];

C = [2 4 3];

D = [0];

Q = fliplr(ctrb(A,B));            % Note: ctrb(A,B) gives Q as [B AB ... A^n-1*B]
                                  % so use fliplr() to get Q as [A^n-1*B ... AB B]

[num, den] = ss2tf(A,B,C,D);

% Create A_bar based on the coefficients of the 'den' of the TF
A_bar = [0 1 0;
         0 0 1;
        -1 -3 -3];

B_bar = [0;
         0;
         1];

% Create C_bar based on the coefficients of the 'num' of the TF
% C_bar = [b0 b1 ... bm 0 ... 0] and has same number of columns as A_bar
C_bar = [5 12 8];
     
Q_bar = fliplr(ctrb(A_bar, B_bar));

P = Q_bar*inv(Q);

ev_desired = [-2 -2.001 -2.0001];

K_bar = place(A_bar, B_bar, ev_desired);

K = K_bar*P;

%% Checking answer

A_check = inv(P)*A_bar*P;
B_check = inv(P)*B_bar;
C_check = C_bar*P;
