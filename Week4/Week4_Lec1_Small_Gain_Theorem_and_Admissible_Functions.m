clear all; 
close all; 

%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

s = tf('s');

G = -0.5/(s+1); % G(s)is proper, stable, and has a magnitude < 1

figure
bode(G);
title('Bode plot of G(s)');

admis1 = 0.5; % constant admissible function

admis2 = 0.5/(s+1); % low pass filter admissible function

T_delay = 1;
admis3 = 0.5*exp(-s*T_delay); % time delay admissible function

fb_admis1 = feedback(G, admis1);
figure
subplot(2,2,1); % Creating a 2x2 grid for multiple plots, and plotting bode in the 1st index
bode(G*admis1); % Note: when doing the bode and nyquist plots, use G(s) * delta(s)
title('G*admis1 bode plot');
subplot(2,2,[2 4]); % Plotting nyquist in the 2nd and 4th index of the subplot grid
nyquist1(G*admis1); % Note: use nyquist1 to plot rational transfer functions
title('G*admis1 nyquist plot');
subplot(2,2,3); % Plotting the step response in the 3rd index of the subplot grid
step(fb_admis1);
title('Step response with admis1');

fb_admis2 = feedback(G, admis2);
figure
subplot(2,2,1);
bode(G*admis2);
title('G*admis2 bode plot');
subplot(2,2,[2 4]);
nyquist1(G*admis2);
title('G*admis2 nyquist plot');
subplot(2,2,3);
step(fb_admis2);
title('Step response with admis2');

fb_admis3 = feedback(G, admis3);
figure
subplot(2,2,1);
bode(G*admis3);
title('G*admis3 bode plot');
subplot(2,2,[2 4]);
nyquist(G*admis3); % Note: use nyquist to plot non-rational tranfer functions since e^-sT is not rational
title('G*admis3 nyquist plot');
subplot(2,2,3);
step(fb_admis3);
title('Step response with admis3');



