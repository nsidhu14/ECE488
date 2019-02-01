clear all;
close all; 

s = tf('s');

%% Uncoupled System (US)
% Kd and Kd gains for critical damping
Kp1 = 3.025;
Kd1 = 1;

Kp2 = 3.025;
Kd2 = 1;

plant_US = [10/((s+1)*s) 0; 0 10/((s+1)*s)];
controller_US = [Kp1+Kd1*s 0; 0 Kp2+Kd2*s];

sys_US = feedback(plant_US*controller_US,eye(2));
figure 
step(sys_US) % this plots the response when r1=step & r2=0 and when r1=0 & r2=step
             % Since this is not a second order system, hence its critical
             % damping step response does not look exactly as a second
             % order critical damped system
title('Uncoupled System Step Response');

%% Coupled System (CS)
plant_CS = [10/((s+1)*s) 1/(s*(s+100)); 1/(s*(s+100)) 10/((s+1)*s)];
controller_CS = [Kp1+Kd1*s 0; 0 Kp2+Kd2*s];

sys_CS = feedback(plant_CS*controller_CS,eye(2));
figure 
step(sys_CS) % this plots the response when r1=step & r2=0 and when r1=0 & r2=step
title('Coupled System Step Response');

