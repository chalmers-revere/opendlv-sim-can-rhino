

H = [1, 0; 0, 1];
f2 = [1.89845, -201.237];
A_n_and = [47.3435, 0.47895; 0, 0; -98.5256, 0];
b_n_and = [13.6186; 1; -23.8717]; 
lb = [0.226381; -1.8];
ub = [0.366381; 1.8];
optoption_1 = optimset('Display', 'off', 'TolFun', 1e-10);
[x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], lb, ub, [], optoption_1)





1.9543
rtA_new 0: 47.3435
rtA_new 1: 0.47895
1.95433
rtA_new 2: 0
rtA_new 3: 0
1.95434
rtA_new 4: -98.5256
rtA_new 5: 0
1.95436
rtb_new 0: 13.6186
1.95438
rtb_new 1: 1
1.95439
rtb_new 2: -23.8717
1.95444
rtlb:  0.226381, -1.8
rtub:  0.366381, 1.8
f2:  1.89845, -201.237
A_n_and.size():  3
 qp2 is solvable? 1
 
 
 
1.95239
rtA_new 0: 47.4927
rtA_new 1: 0.476873
1.95241
rtA_new 2: 0
rtA_new 3: 0
1.95243
rtA_new 4: -98.1211
rtA_new 5: 0
1.95244
rtA_new 6: -63.6516
rtA_new 7: 0.436733
1.95245
rtb_new 0: 13.6952
1.95247
rtb_new 1: 1
1.95249
rtb_new 2: -25.1394
1.9525
rtb_new 3: -18.8128
1.95256
rtlb:  0.226803, -1.8
rtub:  0.366803, 1.8
f2:  2.52367, -266.751
A_n_and.size():  4
 qp2 is solvable? 1
 
 
 