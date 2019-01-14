

H = [1, 0; 0, 1];
f2 = [0; 0];
A_n_and = [10, 21; 0, 1.3];
b_n_and = [-51000; -0]; 
lb = [-100; -10];
ub = [100; 10];
optoption_1 = optimset('Display', 'off', 'TolFun', 1e-10);
[x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], lb, ub, [], optoption_1)