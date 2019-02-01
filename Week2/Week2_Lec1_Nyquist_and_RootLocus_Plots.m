clear all; 
close all; 

s = tf('s');

%% Generating Nyquist Plot
G = 10/(s*(s+2));

figure
nyquist(G);
title('nyquist plot');

figure 
nyquist1(G);
title('nyquist1 plot');

%% Generating Root Locus Plot
n_s = 2;
d_s = s^4 + 3*s^3 + 6*s^2 + s;

sys = n_s/d_s;

figure
rlocus(sys);
title('root locus plot');

% Can increase resolution of the root locus plot by specifying a range of
% 'k' values to plot
k = 0:0.01:17/18;
figure
rlocus(sys,k);
title('root locus plot with k range');


