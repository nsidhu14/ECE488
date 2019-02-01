%% Function to create the linear state-space system

function dxdt = Func_ss_sys_linear(t,x,A,B,u) % pass in the symbolic 'u' and convert it to a number type
   u = double([0; 0]);                        % set tau1 = tau2 = 0 
   dxdt = A*x + B*u;
end