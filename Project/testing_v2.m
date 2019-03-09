clear all; 
close all;

m1 = 0.375;
m2 = 0.375;
l1 = 0.15;
l2 = 0.15;
c1 = 2;
c2 = 2;
g = 3.7;

x_boundary = [0 0 0.22 0.22; 0 0.22 0.22 0]; 
y_boundary = [0 0.22 0.22 0; 0.22 0.22 0 0]; 

lin_space_points = 10;

via_pt_row_curr = 0;
via_pt_col_curr = 0;
via_pt_row_next = 0;
via_pt_col_next = 0;

ev_desired_F = [-1 -1 -1 -1];
ev_A = zeros(4,1,lin_space_points*4);
K = zeros(2,4,lin_space_points*4);
K = complex(K);

% initializing all variables for each via point

A = sym(zeros(4,4,lin_space_points*4));
B = sym(zeros(4,2,lin_space_points*4));
C = sym(zeros(2,4,lin_space_points*4));
D = sym(zeros(2,2,lin_space_points*4));

x_op = zeros(4,1,lin_space_points*4);
tau1_op = zeros(lin_space_points*4,1);
tau2_op = zeros(lin_space_points*4,1);
u_op = zeros(2,1,lin_space_points*4);
x0 = zeros(4,1,lin_space_points*4);

%% IK Model

x_coord = [0.1 0.2 0.2 0.1];
y_coord = [0.2 0.2 0.1 0.1];

q1_angle = zeros(1,4);
q2_angle = zeros(1,4);

for i=1:4
    q1_angle(i) = atan2(y_coord(i),x_coord(i))-acos((x_coord(i)^2 + y_coord(i)^2 + l1^2 - l2^2)/(2*l1*sqrt(x_coord(i)^2+y_coord(i)^2)));
    q2_angle(i) = acos((x_coord(i)^2+y_coord(i)^2-l1^2-l2^2)/(2*l1*l2));
end

%% FK to verify

for i=1:4
    x_FK(i) = l1*cos(q1_angle(i))+l2*cos(q1_angle(i)+q2_angle(i));
    y_FK(i) = l1*sin(q1_angle(i))+l2*sin(q1_angle(i)+q2_angle(i));
end

%% Total angular displacement

TAD = abs(q1_angle(4)-q1_angle(3))+abs(q1_angle(3)-q1_angle(2))+abs(q1_angle(2)-q1_angle(1))+abs(q1_angle(1)-q1_angle(4))+ abs(q2_angle(4)-q2_angle(3))+abs(q2_angle(3)-q2_angle(2))+abs(q2_angle(2)-q2_angle(1))+abs(q2_angle(1)-q2_angle(4))

%% EOAT Position plot

q1_via_point(1,:) = linspace(q1_angle(1),q1_angle(2),lin_space_points);
q1_via_point(2,:) = linspace(q1_angle(2),q1_angle(3),lin_space_points);
q1_via_point(3,:) = linspace(q1_angle(3),q1_angle(4),lin_space_points);
q1_via_point(4,:) = linspace(q1_angle(4),q1_angle(1),lin_space_points);

q2_via_point(1,:) = linspace(q2_angle(1),q2_angle(2),lin_space_points);
q2_via_point(2,:) = linspace(q2_angle(2),q2_angle(3),lin_space_points);
q2_via_point(3,:) = linspace(q2_angle(3),q2_angle(4),lin_space_points);
q2_via_point(4,:) = linspace(q2_angle(4),q2_angle(1),lin_space_points);

x_plot_via_point = [l1*cos(q1_via_point(1,:))+l2*cos(q1_via_point(1,:)+q2_via_point(1,:)),l1*cos(q1_via_point(2,:))+l2*cos(q1_via_point(2,:)+q2_via_point(2,:)),l1*cos(q1_via_point(3,:))+l2*cos(q1_via_point(3,:)+q2_via_point(3,:)),l1*cos(q1_via_point(4,:))+l2*cos(q1_via_point(4,:)+q2_via_point(4,:))];
y_plot_via_point = [l1*sin(q1_via_point(1,:))+l2*sin(q1_via_point(1,:)+q2_via_point(1,:)),l1*sin(q1_via_point(2,:))+l2*sin(q1_via_point(2,:)+q2_via_point(2,:)),l1*sin(q1_via_point(3,:))+l2*sin(q1_via_point(3,:)+q2_via_point(3,:)),l1*sin(q1_via_point(4,:))+l2*sin(q1_via_point(4,:)+q2_via_point(4,:))];

%% Plotting

figure
hold on
scatter(x_plot_via_point,y_plot_via_point, '.b');
plot(x_boundary, y_boundary, '-r');
plot(x_coord, y_coord, 'og');


%% SS

syms q1 q2 q1_d q2_d q1_dd q2_dd tau1 tau2
eqns = [tau1 == ((m1*l1^2)/3 + (m2*l2^2)/12 + m2*(l1^2+(l2^2)/4+l1*l2*cos(q2)))*q1_dd + ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q2_dd - m2*l1*l2*sin(q2)*q1_d*q2_d - (m2*l1*l2*sin(q2))/2*q2_d^2 + ((m1*l1)/2+m2*l1)*g*cos(q1) + (m2*l2)/2*g*cos(q1+q2) + c1*q1_d , tau2 == ((m2*l2^2)/3+(m2*l1*l2)/2*cos(q2))*q1_dd + (m2*l2^2)/3*q2_dd + (m2*l1*l2*sin(q2))/2*q1_d^2 + (m2*l2)/2*g*cos(q1+q2) + c2*q2_d];
dd_vals = solve(eqns, [q1_dd q2_dd]);
tau_vals = solve(eqns, [tau1 tau2]);

x = [q1; q1_d; q2; q2_d];
x_dot = [q1_d; dd_vals.q1_dd; q2_d; dd_vals.q2_dd];
u = [tau1; tau2];
y = [q1; q2];

for i=1:((lin_space_points*4)-1)
    
    A(:,:,i) = jacobian(x_dot,x);
    B(:,:,i) = jacobian(x_dot,u);
    C(:,:,i) = jacobian(y,x);
    D(:,:,i) = jacobian(y,u);
    
    if (i < lin_space_points)
        via_pt_row_next = 1;
    elseif (i >= lin_space_points) && (i < lin_space_points*2)
        via_pt_row_next = 2;
    elseif (i >= lin_space_points*2) && (i < lin_space_points*3)
        via_pt_row_next = 3;
    else
        via_pt_row_next = 4;
    end
    
    if (mod(i+1,lin_space_points) == 0)
        via_pt_col_next = lin_space_points;
        via_pt_col_curr = via_pt_col_next-1;
        via_pt_row_curr = via_pt_row_next;
    elseif (mod(i+1,lin_space_points) == 1)
        via_pt_col_next = mod(i+1,lin_space_points);
        via_pt_col_curr = lin_space_points;
        via_pt_row_curr = via_pt_row_next-1;
    else
        via_pt_col_next = mod(i+1,lin_space_points);
        via_pt_col_curr = via_pt_col_next-1;
        via_pt_row_curr = via_pt_row_next;
    end

    x_op(:,:,i) = [q1_via_point(via_pt_row_next,via_pt_col_next); 0; q2_via_point(via_pt_row_next,via_pt_col_next); 0];
    tau1_op(i) = double(subs(tau_vals.tau1, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(via_pt_row_next,via_pt_col_next) 0 0 q2_via_point(via_pt_row_next,via_pt_col_next) 0 0]));
    tau2_op(i) = double(subs(tau_vals.tau2, [q1 q1_d q1_dd q2 q2_d q2_dd], [q1_via_point(via_pt_row_next,via_pt_col_next) 0 0 q2_via_point(via_pt_row_next,via_pt_col_next) 0 0]));
    u_op(:,:,i) = [tau1_op(i); tau2_op(i)];

    x0(:,:,i) = [q1_via_point(via_pt_row_curr,via_pt_col_curr); 0; q2_via_point(via_pt_row_curr,via_pt_col_curr); 0];

    % Substituing the equilibrium points and converting from syms to numbers type 
    A(:,:,i) = double(subs(A(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));       % Note: use .' to transpose x and u into one long row vector
    B(:,:,i) = double(subs(B(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
    C(:,:,i) = double(subs(C(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
    D(:,:,i) = double(subs(D(:,:,i), [x.' u.'], [x_op(:,:,i).' u_op(:,:,i).']));
    
    x_next = l1*cos(q1_via_point(via_pt_row_next,via_pt_col_next))+l2*cos(q1_via_point(via_pt_row_next,via_pt_col_next)+q2_via_point(via_pt_row_next,via_pt_col_next));
    y_next = l1*sin(q1_via_point(via_pt_row_next,via_pt_col_next))+l2*sin(q1_via_point(via_pt_row_next,via_pt_col_next)+q2_via_point(via_pt_row_next,via_pt_col_next));
    plot(x_next, y_next, 'om');
    
    ev_A(:,:,i) = eig(A(:,:,i));
    
    for j = 1:size(ev_A(:,:,i),1)
        if ev_A(j,:,i) < 0
            ev_desired_F(j) = ev_A(j,:,i);
        end
    end

    K(:,:,i) = place(double(A(:,:,i)), double(B(:,:,i)), ev_desired_F);
end

%% Verifying - Note: should fail very last one since xo(last) is not being used

for j=1:lin_space_points
    if x0(:,:,j).' == [q1_via_point(1,j), 0, q2_via_point(1,j), 0]
        disp('correct');
    else
        disp('incorrect');
    end
end
