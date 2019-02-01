%% Function to create the state-space system

function dxdt = Week2_Tut_Func_ss_sys(t,x,A,B,u) % pass in the symbolic 'u' and convert it to a number type
    if t<1
        u = 1;
    else
        u = 0;
    end
   dxdt = A*x + B*u;
end

