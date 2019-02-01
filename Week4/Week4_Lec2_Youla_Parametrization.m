clear all;
close all;

%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

%% Stabilizing plant for 3 values of r(s), without considering steady-state error to a step input

s = tf('s');

G = 100/((s+1)*(s-2)*(s+4));             % Plant

[A,B,C,D] = ssdata(G); 

p = [-5,-2.5,-1];                        % Picking 3 stable closed-loop poles. 
                                         % The closed-loop poles are the eigenvalues of A-BK and A-FC (for this case).
                                         % Hence, these poles assure that A-BK and A-FC are stable.

K = place(A,B,p);                        % Finding K such that A-BK is stable

F_transpose = place(A',C',p);            % Finding F_transpose such that A-FC is stable
F = F_transpose';                        % Note: Must transpose since A-FC is actually (A^T) - (C^T)*(F^T)

xp = ss(A-F*C, -F, -K, 0);                          % Create the steady-state system for xp
[A_xp, B_xp, C_xp, D_xp] = ssdata(xp);              % Extracting the steady-state matrices of xp
[num_xp, den_xp] = ss2tf(A_xp, B_xp, C_xp, D_xp);   % Getting the numerator and denominator coeffs for xp as a TF 
xp = tf(num_xp,den_xp);                             % Creating the xp TF

np = ss(A-B*K, B, C-D*K, D);
[A_np, B_np, C_np, D_np] = ssdata(np);
[num_np, den_np] = ss2tf(A_np, B_np, C_np, D_np);
np = tf(num_np,den_np);

yp = ss(A-F*C, -B+F*D, -K, 1);
[A_yp, B_yp, C_yp, D_yp] = ssdata(yp);
[num_yp, den_yp] = ss2tf(A_yp, B_yp, C_yp, D_yp);
yp = tf(num_yp,den_yp);

dp = ss(A-B*K, B, -K, 1);
[A_dp, B_dp, C_dp, D_dp] = ssdata(dp);
[num_dp, den_dp] = ss2tf(A_dp, B_dp, C_dp, D_dp);
dp = tf(num_dp,den_dp);

% First choice of r(s)
% Note: Couldn't use syms 'r' since the step function doesn't accept symbolic, 
%       it only takes the variable 's' which is of type 'tf'.
r1 = 0;
C1 = (xp + r1*dp)/(yp - r1*np);
figure
step(feedback(G*C1, 1));
title('Step response with r(s) = 0');

% Second choice of r(s)
r2 = 1;
C2 = (xp + r2*dp)/(yp - r2*np);
figure
step(feedback(G*C2, 1));
title('Step response with r(s) = 1');

% Third choice of r(s)
r3 = 1/(s+1);
C3 = (xp + r3*dp)/(yp - r3*np);
figure
step(feedback(G*C3, 1));
title('Step response with r(s) = 1/(s+1)');

%% Stabilizing plant for 3 values of r(s), with zero steady-state error to a step input

s = tf('s');

G = 100/((s+1)*(s-2)*(s+4));             % Plant
G_prime = G/s;                           % Augmented plant with 1 pure integrator

[A,B,C,D] = ssdata(G_prime); 

p = [-10,-5,-2.5,-1];                    % Picking 4 stable closed-loop poles. 

K = place(A,B,p);                        

F_transpose = place(A',C',p);           
F = F_transpose';                        

xp = ss(A-F*C, -F, -K, 0);                          
[A_xp, B_xp, C_xp, D_xp] = ssdata(xp);              
[num_xp, den_xp] = ss2tf(A_xp, B_xp, C_xp, D_xp);    
xp = tf(num_xp,den_xp);                             

np = ss(A-B*K, B, C-D*K, D);
[A_np, B_np, C_np, D_np] = ssdata(np);
[num_np, den_np] = ss2tf(A_np, B_np, C_np, D_np);
np = tf(num_np,den_np);

yp = ss(A-F*C, -B+F*D, -K, 1);
[A_yp, B_yp, C_yp, D_yp] = ssdata(yp);
[num_yp, den_yp] = ss2tf(A_yp, B_yp, C_yp, D_yp);
yp = tf(num_yp,den_yp);

dp = ss(A-B*K, B, -K, 1);
[A_dp, B_dp, C_dp, D_dp] = ssdata(dp);
[num_dp, den_dp] = ss2tf(A_dp, B_dp, C_dp, D_dp);
dp = tf(num_dp,den_dp);

% First choice of r(s)
r1 = 0;
C1 = (xp + r1*dp)/(yp - r1*np);
C1_prime = C1/s;                   % Shifting the pure integrator from G_prime to C_prime
figure
step(feedback(G*C1_prime, 1));
title('Zero S.S Step response with r(s) = 0');

% Second choice of r(s)
r2 = 1;
C2 = (xp + r2*dp)/(yp - r2*np);
C2_prime = C2/s;
figure
step(feedback(G*C2_prime, 1));
title('Zero S.S Step response with r(s) = 1');

% Third choice of r(s)
r3 = 1/(s+1);
C3 = (xp + r3*dp)/(yp - r3*np);
C3_prime = C3/s;
figure
step(feedback(G_prime*C3, 1));
title('Zero S.S Step response with r(s) = 1/(s+1)');
