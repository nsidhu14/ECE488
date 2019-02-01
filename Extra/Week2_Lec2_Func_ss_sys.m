%% Function to create the state-space system

function dxdt = Week2_Lec2_Func_ss_sys_linear(t,x,A,B,u) % pass in the symbolic 'u' and convert it to a number type
   u = 0.1;
   dxdt = A*x + B*u;
end