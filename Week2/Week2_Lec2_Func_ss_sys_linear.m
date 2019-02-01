%Navdeep Sidhu 20577393
%Ardalan Abolfazli 20571471
%Haiqiao Chen 20569361

%% Function to create the linear state-space system

function dxdt = Week2_Lec2_Func_ss_sys_linear(t,x,A,B,u) % pass in the symbolic 'u' and convert it to a number type
   u = 0.1;
   dxdt = A*x + B*u;
end