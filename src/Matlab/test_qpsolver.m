

H = [1, 0; 0, 1];
f2 = [3.92793e-06, -0.0201236];
A_n_and = [0, 0; 0, 0; -0.125884, -0.999998];
b_n_and = [1; 1.61438; 148.367]; 
lb = [-0.0700017; -3.8];
ub = [ 0.0699983; 3.8];
optoption_1 = optimset('Display', 'off', 'TolFun', 1e-10);
[x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], lb, ub, [], optoption_1)