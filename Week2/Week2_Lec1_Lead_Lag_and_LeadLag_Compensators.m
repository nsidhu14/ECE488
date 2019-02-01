clear all; 
close all; 

%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

s = tf('s');

G = 1/((s+1)*(s+5)*(s+10));
k = 4950;

%% Solving using Lead Compensator
bode(k*G);
[Gm,Pm] = margin(k*G);

Pm_desired = 40 + 5;

% Itteration #1:
theta_max = Pm_desired - Pm;

syms alpha_sym
eqn = sin(degtorad(theta_max)) == (alpha_sym-1)/(alpha_sym+1);
alpha_val = solve(eqn,alpha_sym); % Note: alpha_val is still a symbolic variable so use double(alpha_val)

negative_shift = -10*log10(double(alpha_val)); %REM: USE log10 !!!

Wc = 41.3; % crossover frequency at negative_shift

% Itteration #2:
% Now at Wc, the new phase margin is -68 degrees, so we update theta_max:
theta_max = Pm_desired - -68;

syms alpha_sym
eqn = sin(degtorad(theta_max)) == (alpha_sym-1)/(alpha_sym+1);
alpha_val = solve(eqn,alpha_sym); % Note: alpha_val is still a symbolic variable so use double(alpha_val)

negative_shift = -10*log10(double(alpha_val));

Wc = 30.3; % crossover frequency at negative_shift, has not converged yet

% Itteration #3:
% Now at Wc, the new phase margin is -60 degrees, so we update theta_max:
theta_max = Pm_desired - -60;

syms alpha_sym
eqn = sin(degtorad(theta_max)) == (alpha_sym-1)/(alpha_sym+1);
alpha_val = solve(eqn,alpha_sym); % Note: alpha_val is still a symbolic variable so use double(alpha_val)

negative_shift = -10*log10(double(alpha_val));

Wc = 35.4; % crossover frequency at negative_shift, Hence we have approximately converged

% Once we have converged we can now find 'a' and 'b'
syms a_sym b_sym
eqns = [Wc == sqrt(a_sym*b_sym), double(alpha_val) == b_sym/a_sym];
vars = [a_sym b_sym];
[a_val, b_val] = solve(eqns, vars);

C = k * ((s/double(a_val(2)) + 1)/(s/double(b_val(2)) + 1));

figure
bode(C*G);
title('Lead Compensated Plot');

[Gm_comp_lead,Pm_comp_lead] = margin(C*G);



%% Solving using Lag Compensator
bode(k*G);

% At our desired phase margin of 45 degrees, our phase should be -135
% degrees which occurs at a frequency of 3.99 rad/s. The gain at this
% frequency is 24.8 dB. Hence 20log(alpha)=24.8:

Wc = 3.99;
alpha = 10^(24.8/20);
a = Wc/10;
b = a/alpha; 

C = k * (((s/a)+1)/((s/b)+1));

figure
bode(C*G);
title('Lag Compensated Plot');

[Gm_comp_lag,Pm_comp_lag] = margin(C*G); % almost correct answer, just need to increase Pm_comp_lag by 3 degrees


%% Solving using Lead-Lag Compensator
bode(k*G);