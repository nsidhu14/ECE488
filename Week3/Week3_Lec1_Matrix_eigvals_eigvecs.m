clear all; 
close all; 

%Navdeep Sidhu 20577393
%Haiqiao Chen 20569361

%% Question 1
% Note: matrix A must be of type 'sym' to use the colspace() function
A1 = sym([2 1 0 -1;
         1 0 -1 -2;
        -1 -1 -1 -1]);

RREF_A1 = rref(A1);

rank_A1 = rank(A1);

colspace_A1 = colspace(A1);

nullspace_A1 = null(A1);

%% Question 2

A2 = sym([1 0 1 2 3;
         1 1 0 3 3;
         2 -1 3 3 6]);
     
rank_A2 = rank(A2);

colspace_A2 = colspace(A2);

nullspace_A2 = null(A2);

%% Question 3
  
% No Solution
% since rank[y|A] > rank[A]
y_no_sol = [1;
            1];
A_no_sol = [1 0 0;
            0 0 0];

% Infinite Solutions
% since rank[y|A] == rank[A] (i.e. at least one solution exists) 
% but the augmented matrix [A|y] has a row of all zeros and the number of variables 
% is greater than the number of non-zero rows, so some of the variables can be
% set to any values and still satisfy the equation.
y_inf_sol = [1;
             0];
A_inf_sol = [1 0 0;
             0 0 0];
 
% For a unique solution we require that the matrix A is consistent and that the number of 
% variables is equal to the number of non-zero rows. Hence, a unique solution does not exist 
% since we have 2 equations and 3 variables, so one of the variables can be any value and still satisfy the equation.

%% Question 4

syms epsilon

A4 = [1 0 0; 
     0 3 -4;
     0 1 -1+epsilon];
 
% Note: cannot solve with epsilon = 0 since we cannot disgonalize (not enough eigenvectors) 

% Solving with epsilon = 0.1
A_01 = double(subs(A4, epsilon, 0.1));
[V_01,D_01] = eig(A_01); % diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, 
                         % so that A*V = V*D.

% Solving with epsilon = 0.01
A_001 = double(subs(A4, epsilon, 0.01));
[V_001,D_001] = eig(A_001); 

%% Question 5

A5 = [2 -2 2;
      1 0 2;
      0 0 2];

[V5,D5] = eig(A5);
 