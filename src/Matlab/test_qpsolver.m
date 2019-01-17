

H = [1, 0; 0, 1];
f2 = [0; 0.456982];
A_n_and = [0, 0; 0, 0; -1.57117, 0.999753];
b_n_and = [1.6; 1.6; 75.9421]; 
lb = [-0.5; -4];
ub = [0.5; 4];
optoption_1 = optimset('Display', 'off', 'TolFun', 1e-10);
[x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], lb, ub, [], optoption_1)