clear all;
close all;

s = tf('s');

%% Creating a peak in T(s) using an underdamped system

P = 1/(s*(s+0.1));

Kd = 1;
Kp = 10;
C = Kd*s + Kp;

T = P*C/(1+P*C);
S = 1/(1+P*C);

figure
bodeargs = bodeoptions('cstprefs');
bodeargs.Grid = 'on';
bode(S,T,bodeargs);
title('Bode plot of underdamped response of S(s) and T(s)');

w0 = 3;                       % From the bode plot of S(s), the resonant peak occurs approximately at w = 3 rad/s
sine_input = w0/(s^2+w0^2);
Y_s = T*sine_input;           % Laplace domain output
figure
impulse(Y_s,10);              % Impulse response of Y_s to get time domain output y(t) to a sine input (2nd arg is duration of plot)
title('Impulse response output to a sine input');                              

%% Creating Overdamped Response

P = 1/(s*(s+0.1));

Kd = 10;
Kp = 1;
C = Kd*s + Kp;

T = P*C/(1+P*C);
S = 1/(1+P*C);

figure
bodeargs = bodeoptions('cstprefs');
bodeargs.Grid = 'on';
bode(S,T,bodeargs);
title('Bode plot of overdamped response of S(s) and T(s)');