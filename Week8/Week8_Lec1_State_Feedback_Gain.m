clear all;
close all;

A = [0 1;
    -2 -2];

B = [1;
    0];

ev_desired = [-1.001 -1]; % Note: Must set one ev slightly above -1 since Matlab doesn't allow all ev's to be the same

K = place(A,B,ev_desired);

% Verify final eigenvalues with state feedback gain
ev_final = eig(A-B*K);
