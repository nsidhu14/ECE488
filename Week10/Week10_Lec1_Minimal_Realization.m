clear all;
close all; 

%% Continuing with solving example in lecture notes for section G - adding the observable part (see bottom part)

A = [0 1 0 0 0 0 0;
    -2 -3 0 0 0 0 0;
    0 0 0 1 0 0 0;
    0 0 -2 -3 0 0 0;
    0 0 0 0 0 1 0;
    0 0 0 0 -2 -3 0;
    0 0 0 0 0 0 0];

B = [0 0;
    1 0;
    0 0;
    1 0;
    0 0;
    0 1;
    0 1];

C = [4 2 0 0 4 2 0 ;
    0 0 1 1 0 0 0];

D = [0 0;
    0 0];

Q = ctrb(A,B);

n = size(A,1);
q = rank(Q);

q_lin_indep_cols = [Q(:,(1:4)) Q(:,6)];
    
n_minus_q_lin_indep_cols = [1 0;
                            0 1;
                            0 0;
                            0 0;
                            0 0;
                            0 0;
                            1 0];
                            
T = [q_lin_indep_cols n_minus_q_lin_indep_cols];
T_inv = inv(T);

A_bar = T_inv*A*T;
B_bar = T_inv*B;
C_bar = C*T;
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

% R_bar_1 = obsv(A_bar_11,C_bar_1);
% rank_R_bar_1 = rank(R_bar_1);
% 
% if rank_R_bar_1 == size(A_bar_11,1)
%     disp('R_bar_1 is observable');
% else
%     disp('R_bar_1 is not observable');
% end

%% Taking the controllable part of the system and extracting the minimal realization (i.e. the controllable and observable part)

A_cont = A_bar_11;
B_cont = B_bar_1;
C_cont = C_bar_1;
D_cont = D_bar;

n_ob = size(A_cont,1);

R = obsv(A_cont, C_cont);

r = rank(R);

r_lin_indep_rows = [R((1:2),:)];
n_ob_minus_r_lin_indep_rows = [0 0 0 0 1;
    0 0 0 1 0;
    0 0 1 0 0];

T_ob = [r_lin_indep_rows;
    n_ob_minus_r_lin_indep_rows];

T_ob_inv = inv(T_ob);

A_bar_ob = T_ob*A_cont*T_ob_inv;
B_bar_ob = T_ob*B_cont;
C_bar_ob = C_cont*T_ob_inv;
D_bar_ob = D_cont;

A_bar_11_ob = A_bar_ob((1:r),(1:r));
B_bar_1_ob = B_bar_ob((1:r),:);
C_bar_1_ob = C_bar_ob(:,(1:r));
D_bar_1 = D_bar_ob;

if (rank(ctrb(A_bar_11_ob,B_bar_1_ob)) == size(A_bar_11_ob,1)) && (rank(obsv(A_bar_11_ob,C_bar_1_ob)) == size(A_bar_11_ob,1))
    disp('The system is minimal')
else
    disp('The system is not yet minimal');
end

