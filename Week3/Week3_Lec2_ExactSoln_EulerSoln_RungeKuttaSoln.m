clear all;
close all;


A = [-2 1;
      0 -2];
  
B = [0;
     1];

C = [1 1];

D = [0];
 
x0 = [0;
      0];

u = 1; % input is a step

%% Finding exact solution

syms s t tau

eAt = ilaplace((s*eye(2)-A)^-1);

f = (subs(eAt, t, tau))^-1;                % eAtau^-1
fInt = int(f, tau, [0 t]);                 % integral term of exact solution

x_exact = eAt*x0 + eAt*fInt*B;

y_exact = C*x_exact + D*u;

figure
fplot(y_exact,[0 10]);                     % use fplot since y_exact contains x_exact which has symbolic var 't'
title('Exact Solution');

%% Finding Euler method solution 

% Initial conditions and setup
delta_T = 0.1;                                % step size
tspan = 0:delta_T:10;                         % the range of t
y = zeros(size(tspan));                       % allocate the result y
x(:,1) = double(subs(x_exact, t, 0));         % the initial x value
y(1) = C*x(:,1) + D*u;                        % the initial y value
n = numel(y);                                 % the number of y values
% The loop to solve the DE
for i=1:n-1
    x_dot(:,i) = A*x(:,i) + B*u;              % x_dot at current time step
    x(:,i+1) = delta_T * x_dot(:,i) + x(:,i); % x at next time step
    y(i+1) = C*x(:,i+1) + D*u;                % y at next time step 
end

figure
plot(tspan, y);
title('Euler Method Solution');

%% Finding Runge-Kutta solution - this is for a constant step input 'u'

delta_T = 0.1;                                     % step size
tspan = 0:delta_T:10;                              % time span
y = zeros(size(tspan)); 
x(:,1) = double(subs(x_exact, t, 0));              % the initial x value
y(1) = C*x(:,1) + D*u;                             % the initial y value
n = numel(y);                                      % the number of y values

Func_x_dot = @(t,x) A*x + B*u_input(t);            % function for x_dot at current itteration,
                                                   % u_input(t) is a function for the input defined at the end of the script                                                   

for i=1:n-1                                        % calculation loop 
    k_1 = Func_x_dot(tspan(i), x(:,i));
    k_2 = Func_x_dot(tspan(i)+0.5*delta_T, x(:,i)+0.5*k_1); % Note: 0.5*k_1 is added to each element in the matrix
    k_3 = Func_x_dot(tspan(i)+0.5*delta_T, x(:,i)+0.5*k_2); % Note: 0.5*k_2 is added to each element in the matrix
    k_4 = Func_x_dot(tspan(i)+delta_T, x(:,i)+k_3);         % Note: k_3 is added to each element in the matrix

    x(:,i+1) = x(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*delta_T;  % x at next time step
    y(i+1) = C*x(:,i+1) + D*u;                                % y at next time step
end

figure
plot(tspan, y);
title('Runge-Kutta Method Solution');

% Function definition for input 'u'
function u_t = u_input(t)
    if t >= 0
        u_t = 1;
    else
        u_t = 0;
    end
end