
%% Solving the mass-spring-damper s.s. eqns using ode45

function [] = Week2_Tut_Func_call_dstate() % Syntax: function [out1,out2, ..., outN] = myfun(in1,in2,in3, ..., inN)
                                           % Note: must save the file name as the function name
          
    tspan = [0 20]; % time span of simulation
    x0 = [0;0]; % initial conditions for x1 and x2
    
    % Reference: [t,y] = ode45(odefun,tspan,y0), where tspan = [t0 tf], integrates the system of differential equations y'=f(t,y) 
    % from t0 to tf with initial conditions y0. Each row in the solution array y corresponds to a value returned 
    % in column vector t.
    
    [t,x] = ode45(@dstate,tspan,x0); %@ is a function handler to pass one function into another
    
    plot(t,x(:,1)); % plotting the results 
    
    function dxdt = dstate(t,x) % this function must accept two inputs (t,x) even if one of the inputs is not used.
        % defining constants
        M = 1;
        B = 1;
        K = 1;
        
        if t<1
            u = 1;
        else
            u = 0;
        end
        
        x1=x(1);
        x2=x(2);
        dxdt=[x2; -B/M*x2 - K/M*x1 + 1/M*u];
    end
end