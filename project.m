clc

% To design a doublet achromat with minimized focal length(f) 
% under the constraints of certain radii of both lens, with coma (c) and spherical aberration (s) corrected.

% Radii constraints
% R1 = -R2, R2 = R3, R1 < 400mm, R4 < -1600mm

% Lens matrials table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_c, n_._1 stands for BK7 and n_._2 stands for F2
n_c_1 = 1.514322;
n_c_2 = 1.615032;
% n_d
n_d_1 = 1.516800;
n_d_2 = 1.620040;
% n_e
n_e_1 = 1.518722;
n_e_2 = 1.624080;
% n_f
n_f_1 = 1.522376;
n_f_2 = 1.632081;
% Abbe number V_e
V_1 = 64.4;
V_2 = 36.6;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms q_1 q_2 f

% step 1: depict focal length with f
f_1 = (V_1 - V_2)*f/V_1;
f_2 = (V_2 - V_1)*f/V_2;

% step 2: put f1 and f2 into c1 and c2, with p1=-1 and p2=2.52
% TBD: explain how come p1=-1 and p2=2.52
p_1 = -1;
p_2 = 2.52;

c_1 = ((2*n_e_1 + 1)*p_1 + ((n_e_1 + 1)/(n_e_1 - 1))*q_1)/(4*(n_e_1)*(f_1^2));
c_2 = ((2*n_e_2 + 1)*p_2 + ((n_e_2 + 1)/(n_e_2 - 1))*q_2)/(4*(n_e_2)*(f_2^2));

% step 3: put f1 and f2 into s1 and s2, with p1=-1 and p2=2.52.
S_1 = -(n_e_1^3 + (n_e_1 + 2)*q_1^2 + (3*n_e_1 + 2)*(n_e_1 - 1)^2*p_1^2 + 4*(n_e_1^2 - 1)*p_1*q_1)/(32*n_e_1*(n_e_1 - 1)^2*f_1^3);
S_2 = -(n_e_2^3 + (n_e_2 + 2)*q_2^2 + (3*n_e_2 + 2)*(n_e_2 - 1)^2*p_2^2 + 4*(n_e_2^2 - 1)*p_2*q_2)/(32*n_e_2*(n_e_2 - 1)^2*f_2^3);

% step 4: derive f=F(q_1, q_2) from c_1 + c_2 + S_1 + S_2 = 0
temp = c_1 + c_2 + S_1 + S_2;
F = solve(temp, f);
