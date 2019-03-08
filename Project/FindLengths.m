%% Finding optimal lengths

%Constraints
%Total mass = 0.75kg
%Material choices(kg/m): 0.85, 1.4, 2.5, 7

%Logic is to minimize total angular displacement - minimizes time and
%energy (Don't care enough to account for friction and controls stuff)

%Find optimal lengths based on different material configs, then compare TAD
%and total energy

%Find range of possible lengths to reach all points and meet mass
%constraint

%Case 1: both lightest

%0.75 = 0.85L1 + 0.85L2

%% Lengths

clear all;
close all;

% L1 = 10*sqrt(0.1^2*2);
L1 = 0.5;
L2= (0.75-0.85*L1)/0.85;

m = 0.85*(L1+L2)

%% IK Model

x = [0.2 0.2 0.1 0.1];
y = [0.2 0.1 0.1 0.2];

Q1=zeros(1,4);
Q2=zeros(1,4);
for i=1:4
    Q1(i) = atan2(y(i),x(i))-acos((x(i)^2 + y(i)^2 + L1^2 - L2^2)/(2*L1*sqrt(x(i)^2+y(i)^2)));
    Q2(i) = acos((x(i)^2+y(i)^2-L1^2-L2^2)/(2*L1*L2));
end
    
%% FK to verify

for i=1:4
    x_FK(i) = L1*cos(Q1(i))+L2*cos(Q1(i)+Q2(i));
    y_FK(i) = L1*sin(Q1(i))+L2*sin(Q1(i)+Q2(i));
end

%% Total angular displacement

TAD = abs(Q1(4)-Q1(3))+abs(Q1(3)-Q1(2))+abs(Q1(2)-Q1(1))+abs(Q1(1)-Q1(4))+ abs(Q2(4)-Q2(3))+abs(Q2(3)-Q2(2))+abs(Q2(2)-Q2(1))+abs(Q2(1)-Q2(4))

%% EOAT Position plot
Q1_plot(1,:) = linspace(Q1(1),Q1(2));
Q1_plot(2,:) = linspace(Q1(2),Q1(3));
Q1_plot(3,:) = linspace(Q1(3),Q1(4));
Q1_plot(4,:) = linspace(Q1(4),Q1(1));

Q2_plot(1,:) = linspace(Q2(1),Q2(2));
Q2_plot(2,:) = linspace(Q2(2),Q2(3));
Q2_plot(3,:) = linspace(Q2(3),Q2(4));
Q2_plot(4,:) = linspace(Q2(4),Q2(1));

x_plot = [L1*cos(Q1_plot(1,:))+L2*cos(Q1_plot(1,:)+Q2_plot(1,:)),L1*cos(Q1_plot(2,:))+L2*cos(Q1_plot(2,:)+Q2_plot(2,:)),L1*cos(Q1_plot(3,:))+L2*cos(Q1_plot(3,:)+Q2_plot(3,:)),L1*cos(Q1_plot(4,:))+L2*cos(Q1_plot(4,:)+Q2_plot(4,:))];
y_plot = [L1*sin(Q1_plot(1,:))+L2*sin(Q1_plot(1,:)+Q2_plot(1,:)),L1*sin(Q1_plot(2,:))+L2*sin(Q1_plot(2,:)+Q2_plot(2,:)),L1*sin(Q1_plot(3,:))+L2*sin(Q1_plot(3,:)+Q2_plot(3,:)),L1*sin(Q1_plot(4,:))+L2*sin(Q1_plot(4,:)+Q2_plot(4,:))];

scatter(x_plot,y_plot);