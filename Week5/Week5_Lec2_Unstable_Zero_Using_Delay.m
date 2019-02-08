clear all;
close all;

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361

%% Critically damped PD controler without delay

s = tf('s');

P = 1/(s*(s+10));

Kd = 1;
Kp = 30.25;
C = Kd*s + Kp;

T = P*C/(1+P*C);
S = 1/(1+P*C);

figure
step(T);
title('Step response of closed-loop system without delay');

figure
hold on
bodeargs = bodeoptions('cstprefs');
bodeargs.Grid = 'on';
bode(S,T,bodeargs);
title('Bode plot of S(s) and T(s) without delay');
hold off

%% Critically damped PD controler with introducing various delays

% Delay
Td1 = 0.1;
Td2 = 0.05;
Td3 = 0.01;

% First order Pade approximation of the delay
delay1 = (1-(Td1*s)/2)/(1+(Td1*s)/2);
delay2 = (1-(Td2*s)/2)/(1+(Td2*s)/2);
delay3 = (1-(Td3*s)/2)/(1+(Td3*s)/2);

% Complimentary sensitivity function (i.e. closed-loop transfer function)
T1 = P*C*delay1/(1+P*C*delay1);
T2 = P*C*delay2/(1+P*C*delay2);
T3 = P*C*delay3/(1+P*C*delay3);

% Sensitivity function
S1 = 1/(1+P*C*delay1);
S2 = 1/(1+P*C*delay2);
S3 = 1/(1+P*C*delay3);

figure
hold on
step(T1);
step(T2);
step(T3);
title('Step response of closed-loop system with various delays');
legend('Td=0.1','Td=0.05','Td=0.01');
hold off

figure
hold on
bodeargs = bodeoptions('cstprefs');
bodeargs.Grid = 'on';
bode(S1,T1,bodeargs);
bode(S2,T2,bodeargs);
bode(S3,T3,bodeargs);
title('Bode plot of S(s) and T(s) with various delays');
legend('S: Td=0.1','T: Td=0.1','S: Td=0.05','T: Td=0.05','S: Td=0.01','T: Td=0.01');
hold off