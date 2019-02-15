clear all;
close all;

%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

%% Plotting closed-loop step response w/ diagonalized plant (3 cases)

s = tf('s');

P = [1/(s*(s+1)) 1/(s*(s+10));
     1/(s*(s+10)) 1/(s*(s+1))];
 
P_diag = [1/(s*(s+1)) 0;
          0 1/(s*(s+1))];

% All closed-loop poles at s = -0.5
Kp1 = 0.25;
Kd1 = 0; 
C1 = [Kp1+Kd1*s/(0.1*s+1) 0;
      0 Kp1+Kd1*s/(0.1*s+1)];
  
sys1_diag = feedback(P_diag*C1, eye(2));

figure
step(sys1_diag);
title('Closed-loop step response w/ diagonalized plant (All poles at s=-0.5)'); 

% All closed-loop poles at s = -1
Kp2 = 1;
Kd2 = 1;
C2 = [Kp2+Kd2*s/(0.1*s+1) 0;
      0 Kp2+Kd2*s/(0.1*s+1)];
  
sys2_diag = feedback(P_diag*C2, eye(2));

figure
step(sys2_diag);
title('Closed-loop step response w/ diagonalized plant (All poles at s=-1)'); 

% All closed-loop poles at s = -5
Kp3 = 25;
Kd3 = 9; 
C3 = [Kp3+Kd3*s/(0.1*s+1) 0;
      0 Kp3+Kd3*s/(0.1*s+1)];
  
sys3_diag = feedback(P_diag*C3, eye(2));

figure
step(sys3_diag);
title('Closed-loop step response w/ diagonalized plant (All poles at s=-5)'); 

%% Plotting closed-loop step response w/ actual plant (3 cases)

% All closed-loop poles at s = -0.5  
sys1_actual = feedback(P*C1, eye(2));

figure
step(sys1_actual);
title('Closed-loop step response w/ actual plant (All poles at s=-0.5)'); 

% All closed-loop poles at s = -1  
sys2_actual = feedback(P*C2, eye(2));

figure
step(sys2_actual);
title('Closed-loop step response w/ actual plant (All poles at s=-1)'); 

% All closed-loop poles at s = -5  
sys3_actual = feedback(P*C3, eye(2));

figure
step(sys3_actual);
title('Closed-loop step response w/ actual plant (All poles at s=-5)'); 

%% Decoupling the plant using T

P_low_freq = [1 1/10;
              1/10 1];

T = inv(P_low_freq);

C1_decoupled = T*C1;
C2_decoupled = T*C2;
C3_decoupled = T*C3;

sys1_decoupled = feedback(P*C1_decoupled, eye(2));
sys2_decoupled = feedback(P*C2_decoupled, eye(2));
sys3_decoupled = feedback(P*C3_decoupled, eye(2));

figure
step(sys1_decoupled);
title('Closed-loop step response w/ decoupling T (All poles at s=-0.5)');

figure
step(sys2_decoupled);
title('Closed-loop step response w/ decoupling T (All poles at s=-1)');

figure
step(sys3_decoupled);
title('Closed-loop step response w/ decoupling T (All poles at s=-5)');