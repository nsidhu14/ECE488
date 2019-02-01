A = [ 0 1 0
    980 0 -2.8
    0 0 -100];

B = [0
    0
    100];

C = [1 0 0]; 

%Stability
poles = eig(A);

t = 0:0.01:2;
u = zeros(size(t));
x0 = [0.01 0 0];

sys = ss(A,B,C,0);

% [y,t,x] = lsim(sys,u,t,x0);
% plot(t,y)
% title('Open-Loop Response to Non-Zero Initial Condition')
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')


%Controllability and Observability
Wc = ctrb(A,B);
rank_Wc = rank(Wc);

Wo = obsv(A,C);
rank_Wo = rank(Wo);

%Control Design Using Pole Placement
p1 = -20+20i;
p2 = -20-20i;
p3 = -100; 

K = place(A,B,[p1 p2 p3]);
Nbar = rscale(sys,K)
sys_cl = ss(A-B*K,B,C,0);

% [y,t,x] = lsim(sys_cl,u,t,x0);
% plot(t,y)
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')

%Introducing the Reference Input
t = 0:0.01:2;
u = 0.001*ones(size(t));

sys_cl = ss(A-B*K,B,C,0);

% [y,t,x] = lsim(sys_cl,u,t);
% plot(t,y)
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')
% axis([0 2 -4E-6 0])

% [y,t,x] = lsim(sys_cl,Nbar*u,t)
% plot(t,y)
% title('Linear Simulation Results (with Nbar)')
% xlabel('Time (sec)')
% ylabel('Ball Position (m)')
% axis([0 2 0 1.2*10^-3])

%Observer Design
op1 = -100;
op2 = -101;
op3 = -102;

L = place(A',C',[op1 op2 op3])';

At = [ A-B*K             B*K
       zeros(size(A))    A-L*C ];

Bt = [    B*Nbar
       zeros(size(B)) ];

Ct = [ C    zeros(size(C)) ];

sys = ss(At,Bt,Ct,0);
[y,t,x] = lsim(sys,zeros(size(t)),t,[x0 x0]);
plot(t,y)
title('Linear Simulation Results (with observer)')
xlabel('Time (sec)')
ylabel('Ball Position (m)')


t = 0:1E-6:0.1;
x0 = [0.01 0.5 -5];
[y,t,x] = lsim(sys,zeros(size(t)),t,[x0 x0]);

n = 3;
e = x(:,n+1:end);
x = x(:,1:n);
x_est = x - e;

% Save state variables explicitly to aid in plotting
h = x(:,1); h_dot = x(:,2); i = x(:,3);
h_est = x_est(:,1); h_dot_est = x_est(:,2); i_est = x_est(:,3);

plot(t,h,'-r',t,h_est,':r',t,h_dot,'-b',t,h_dot_est,':b',t,i,'-g',t,i_est,':g')
legend('h','h_{est}','hdot','hdot_{est}','i','i_{est}')
xlabel('Time (sec)')