syms a b
freq_crossover = 5;
alpha = 5.8;
eqns = [freq_crossover == sqrt(a*b), alpha == b/a];
vars = [a b];
[a_val, b_val] = solve(eqns, vars);
vpa(a_val)
vpa(b_val)