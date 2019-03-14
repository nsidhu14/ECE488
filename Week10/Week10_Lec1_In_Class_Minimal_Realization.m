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

Q = ctrb(A,B);

q = rank(Q);
n_cont = size(A,1);

lin_indep_cols_of_Q = orth(Q);  % Note: orth(Q) returns linearly independant columns which span Q

q_lin_indep_cols = lin_indep_cols_of_Q(:,(1:q));    % we requrie q lin. indep. columns from orth(Q)

n_minus_q_lin_indep_cols = [1 0;
    0 1;
    0 0;
    0 1;
    0 0;
    1 0;
    0 0;
    0 0];

T_cont = [q_lin_indep_cols n_minus_q_lin_indep_cols];

T_cont_inv = inv(T_cont);

A_bar = T_cont_inv*A*T_cont;
B_bar = T_cont_inv*B;
C_bar = C*T_cont;
D_bar = D;

A_bar_11 = A_bar((1:q),(1:q));
B_bar_1 = B_bar((1:q),:);
C_bar_1 = C_bar(:,(1:q));
D_bar_1 = D;

Q_bar_1 = ctrb(A_bar_11,B_bar_1);
rank_Q_bar_1 = rank(Q_bar_1);

if rank_Q_bar_1 == size(A_bar_11,1)
    disp('Q_bar_1 is controllable');
else
    disp('Q_bar_1 is not controllable');
end

R_bar_1 = obsv(A_bar_11,C_bar_1);
rank_R_bar_1 = rank(R_bar_1);

if rank_R_bar_1 == size(A_bar_11,1)
    disp('R_bar_1 is observable');
else
    disp('R_bar_1 is not observable');
end

%% Checking if the final system is equivilant to the original system

epsilon = 10^-8;    % using epsilion to account for very small numbers

if (abs(D - D_bar) < (epsilon + zeros(size(D))))
    D_check = 1;
else
    D_check = 0;
end

if (abs(C*B - C_bar_1*B_bar_1) < (epsilon + zeros(size(C*B))))
    CB_check = 1;
else
    CB_check = 0;
end

if (abs(C*A*B - C_bar_1*A_bar_11*B_bar_1) < (epsilon + zeros(size(C*A*B))))
    CAB_check = 1;
else
    CAB_check = 0;
end

if (abs(C*(A)^2*B - C_bar_1*(A_bar_11)^2*B_bar_1) < (epsilon + zeros(size(C*(A)^2*B))))
    CA2B_check = 1;
else
    CA2B_check = 0;
end

if (abs(C*(A)^3*B - C_bar_1*(A_bar_11)^3*B_bar_1) < (epsilon + zeros(size(C*(A)^3*B))))
    CA3B_check = 1;
else
    CA3B_check = 0;
end

if (abs(C*(A)^4*B - C_bar_1*(A_bar_11)^4*B_bar_1) < (epsilon + zeros(size(C*(A)^4*B))))
    CA4B_check = 1;
else
    CA4B_check = 0;
end

if (abs(C*(A)^5*B - C_bar_1*(A_bar_11)^5*B_bar_1) < (epsilon + zeros(size(C*(A)^5*B))))
    CA5B_check = 1;
else
    CA5B_check = 0;
end

if ((D_check == 1) && (CB_check == 1) && (CAB_check == 1) && (CA2B_check == 1) && (CA3B_check == 1) && (CA4B_check == 1) && (CA5B_check == 1))
    disp('The two transfer functions are the same');
else
    disp('The two transfer functions are not the same');
end
