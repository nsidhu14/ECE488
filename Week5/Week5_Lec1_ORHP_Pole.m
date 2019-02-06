clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361

s = tf('s');

%% Case 1: b = 0.1

b = 0.1;

P = 1/(s*(s-b));

a = 3+b;
Kd = 3+a*b;
Kp = 1;

C = (Kd*s+Kp)/(s+a);

T1 = P*C/(1+P*C);
S1 = 1/(1+P*C);

figure
step(T1);
title('Step response with b=0.1');
figure 
bode(S1)
title('Bode plot of S(s) with b=0.1');

%% Case 2: b = 0.5

b = 0.5;

P = 1/(s*(s-b));

a = 3+b;
Kd = 3+a*b;
Kp = 1;

C = (Kd*s+Kp)/(s+a);

T2 = P*C/(1+P*C);
S2 = 1/(1+P*C);

figure
step(T2);
title('Step response with b=0.5');
figure 
bode(S2)
title('Bode plot of S(s) with b=0.5');

%% Case 3: b = 1

b = 1;

P = 1/(s*(s-b));

a = 3+b;
Kd = 3+a*b;
Kp = 1;

C = (Kd*s+Kp)/(s+a);

T3 = P*C/(1+P*C);
S3 = 1/(1+P*C);

figure
step(T3);
title('Step response with b=1');
figure 
bode(S3)
title('Bode plot of S(s) with b=1');

%% Plotting all cases at the same time

figure 
hold on 
step(T1);
step(T2);
step(T3);
title('Step response for all three cases');
legend('b=0.1','b=0.5','b=1');

figure 
hold on 
bode(S1);
bode(S2);
bode(S3);
title('Bode plot for all three cases');
legend('b=0.1','b=0.5','b=1');