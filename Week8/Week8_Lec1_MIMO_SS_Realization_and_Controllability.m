clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361
%Ardalan Abolfazli 20571471

%% Finding SS realization of MIMO system

s = tf('s');

G = [1/(s*(s+1)) 1/(s*(s+10));
    1/(s*(s+10)) 1/(s*(s+1))];

[num11,den11] = tfdata(G(1,1)); % Note: num11 and den11 are cell arrays
[num21,den21] = tfdata(G(2,1)); 
[num12,den12] = tfdata(G(1,2)); 
[num22,den22] = tfdata(G(2,2)); 

num11 = cell2mat(num11);        % Must use cell2mat() to convert the cell array to an ordinary matrix
num21 = cell2mat(num21);
num12 = cell2mat(num12);
num22 = cell2mat(num22);
den11 = cell2mat(den11);
den21 = cell2mat(den21);
den12 = cell2mat(den12);
den22 = cell2mat(den22);

% SS realization of each element
[A11,B11,C11,D11] = tf2ss(num11, den11);
[A21,B21,C21,D21] = tf2ss(num21, den21);
[A12,B12,C12,D12] = tf2ss(num12, den12);
[A22,B22,C22,D22] = tf2ss(num22, den22);

% SS realization for the entire system
A = [A11 zeros(2) zeros(2) zeros(2);
    zeros(2) A21 zeros(2) zeros(2);
    zeros(2) zeros(2) A12 zeros(2);
    zeros(2) zeros(2) zeros(2) A22];

B = [B11 zeros(2,1);
    B21 zeros(2,1);
    zeros(2,1) B12;
    zeros(2,1) B22];

C = [C11 zeros(1,2) C12 zeros(1,2);
    zeros(1,2) C21 zeros(1,2) C22];

D = [D11 D12;
    D21 D22];

%% Checking answer by re-creating the TF using A,B,C,D and seeing if its approx. G(s) 

G_check = C*inv(s*eye(size(A))-A)*B+D;
G_check_mr = minreal(G_check); % minreal() cancels pole-zero pairs in the TF

% Since G_check_mr basically looks the same as G(s), thus the SS
% realization was correct.

%% Checking Controllability 

% since A is 8x8, thus Q has max of A^7*B
Q = [A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B A*B B];
rank_Q = rank(Q);

% Or use ctrb command
Q_ctrb = ctrb(A,B);
rank_Q_ctrb = rank(Q_ctrb);

% Since Q is 8x16 and the rank of Q is 6, hence Q is not full rank (i.e.
% does not have a rank of 8). Thus, the system is not controllable.  



