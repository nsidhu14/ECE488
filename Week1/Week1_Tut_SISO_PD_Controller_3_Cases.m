clear all;
close all;

% Define the laplace variable 's'
s = tf('s');

% Plant
P = 10/(s^2+10*s);

%Bode plot of plant
figure
bode(P)
title('Plant Bode Plot');

% Gain and phase margin of plant
[Gm,Pm] = margin(P)

%-------------------------------------------------------------------------
% Pd controller C = Kp + Kd*s
% char poly = s^2 + (10+10*Kd)*s + 10*Kp
% 1) Critically damped PD controller:
% pick Kp and Kd such that char poly has b^2 = 4*a*c
% b^2 = 100 + 200*Kd + 100*Kd^2
% 4*a*c = 40*Kp

% By setting Kp = 10, we can solve for Kd = 1 and hence achieve a
% critically damped system

Kp = 10;
Kd = 1;

C = Kp + Kd*s;

% Closed loop system, REM do P*C in the feedback
sys_crit = feedback(P*C,1);

% Critically damped bode plot
figure 
bode(sys_crit)
title('Critically Damped Bode Plot');

% Critically damped step response
figure
step(sys_crit)
title('Critically Damped Step Response');

%-------------------------------------------------------------------------
% 2) Underdamped PD controller:
% pick Kp and Kd such that char poly has b^2 < 4*a*c
% b^2 = 100 + 200*Kd + 100*Kd^2
% 4*a*c = 40*Kp

% By setting Kp = 10, we can solve for Kd < 1 and hence achieve an
% underdamped system

Kp = 10;
Kd = 0.01;

C = Kp + Kd*s;

% Closed loop system
sys_under = feedback(P*C,1);

% Underdamped bode plot
figure 
bode(sys_under)
title('Underdamped Bode Plot');

% Underdamped step response
figure 
step(sys_under)
title('Underdamped Step Response');

%-------------------------------------------------------------------------
% 3) Overdamped PD controller:
% pick Kp and Kd such that char poly has b^2 > 4*a*c
% b^2 = 100 + 200*Kd + 100*Kd^2
% 4*a*c = 40*Kp

% By setting Kp = 10, we can solve for Kd > 1 and hence achieve an
% underdamped system

Kp = 10;
Kd = 100;

C = Kp + Kd*s;

% Closed loop system
sys_over = feedback(P*C,1);

% Overdamped bode plot
figure 
bode(sys_over)
title('Overdamped Bode Plot');

% Overdamped step response
figure 
step(sys_over)
title('Overdamped Step Response');

%-------------------------------------------------------------------------
% Plotting all 3 case's step responses on one plot
figure 
hold on
step(sys_crit)
step(sys_under)
step(sys_over)
title('Step Response of All Three Cases');
hold off

%------------------------------------------------------------------------
% In reality we can't have controller C = Kp + Kd*s since we don't have
% infinite bandwidth. Hence, we must augment the controller with a low pass
% filter alpha/(s+alpha)(i.e. 1/((s/alpha)+1) so that at low freq where s = j0 
% the effect is approx. negligiable, and at high freq where s = j*infinity we
% get a cut-off freq at alpha. Now when we set alpha sufficiently high such
% as 1000, we see that our systems behave just as they did w/o the low pass
% filter.

alpha = 1000;
low_pass_filter = 1/((s/alpha)+1);

% Critical damping with low pass filter
Kp = 10;
Kd = 1;
C = Kp + Kd*s;
sys_crit = feedback(P*low_pass_filter*C,1);

% Underdamping with low pass filter
Kp = 10;
Kd = 0.01;
C = Kp + Kd*s;
sys_under = feedback(P*low_pass_filter*C,1);

% Overdamping with low pass filter
Kp = 10;
Kd = 100;
C = Kp + Kd*s;
sys_over = feedback(P*low_pass_filter*C,1);

% Plotting all step responses with low pass filter
figure 
hold on
step(sys_crit)
step(sys_under)
step(sys_over)
title('Step Response of All Three Cases With Low Pass Filter');
hold off


