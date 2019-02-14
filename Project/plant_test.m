clear all
close all

c1=1;
c2=1;
g=3.7;
m1=1;
m2=1;
l1=1;
l2=1;

x01=[0 0 0 0]';
yout1=[0 0];
for t=0:0.001:2
[tout1,xout1]=ode45(@(t,x)simulatorofficial(t,x,[1;1],l1,l2,m1,m2,g,c1,c2),[t t+0.001],x01); %put your plant here
y1=[1 0 0 0;0 0 1 0]*xout1(end,:)';
yout1=[yout1;y1'];
x01=xout1(end,:)';
t
end
tspan=0:0.001:2;
yout1=yout1(2:end,:);
plot(tspan,yout1);
legend('q1','q2')

x02=[pi/2 0 0 0]';
yout2=[0 0];
for t=0:0.001:2
[tout2,xout2]=ode45(@(t,x)simulatorofficial(t,x,[sin(t);cos(t)],l1,l2,m1,m2,g,c1,c2),[t t+0.001],x02); %put your plant here
y2=[1 0 0 0;0 0 1 0]*xout2(end,:)';
yout2=[yout2;y2'];
x02=xout2(end,:)';
t
end
tspan=0:0.001:2;
figure
yout2=yout2(2:end,:);
plot(tspan,yout2);
legend('q1','q2')