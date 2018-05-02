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

% step 5: add penalty
% R1 < 400mm, R4 < -1600mm
disp('Start to do penalty');
R_1 = 2*F*((n_e_1 - 1) - (n_f_1 - n_c_1)*(n_e_2 - 1)/(n_f_2 - n_c_2))/(q_1 + 1);
R_3 = -((n_f_2 - n_c_2)*(q_1 + 1)*R_1)/((n_f_1 - n_c_1)*(q_2 + 1));
R_4 = (q_2 + 1)*R_3/(q_2 - 1);

Penal_1 = R_1 - 400;
Penal_2 = R_4 + 1600;
R = 1;

F = F + R*heaviside(Penal_1)*(Penal_1^2) + R*heaviside(Penal_2)*(Penal_2^2);

% step 6: use Conjugate Gradient Method to minimize F
g = gradient(F, [q_1, q_2]);
fcontour(F, 'LevelStep', 5);
grid on;
axis([-100 200 -100 200]);
axis square
hold on

F_fun = matlabFunction(F);

% derivative
s1 = matlabFunction(g(1));
s2 = matlabFunction(g(2));

% initial value for CG
q_1_cur = 400;
q_2_cur = 400;
F_cur = F_fun(q_1_cur, q_2_cur);
q_1_pre = Inf;
q_2_pre = Inf;
i = 1;

disp(['q_1_cur: ' num2str(q_1_cur)]);
disp(['q_2_cur: ' num2str(q_2_cur)]);
disp(['F_cur: ' num2str(F_cur)]);
scatter(q_1_cur, q_2_cur, 'r');
text(q_1_cur, q_2_cur, 'Initial point');

while ((abs((q_1_cur - q_1_pre)/q_1_cur) > 0.01) || (abs((q_2_cur - q_2_pre)/q_2_cur) > 0.01))
    syms d

    if(i == 1)
        q_1_d = q_1_cur + d*(-1)*s1(q_1_cur, q_2_cur);
        q_2_d = q_2_cur + d*(-1)*s2(q_1_cur, q_2_cur);

    else
        tmp1 = s1(q_1_cur, q_2_cur)^2 + s2(q_1_cur, q_2_cur)^2;
        tmp2 = s1(q_1_pre, q_2_pre)^2 + s2(q_1_pre, q_2_pre)^2;
        q_1_d = q_1_cur + d*((-1)*s1(q_1_cur, q_2_cur) + (tmp1/tmp2)*s1(q_1_cur, q_2_cur));
        q_2_d = q_1_cur + d*((-1)*s2(q_1_cur, q_2_cur) + (tmp1/tmp2)*s2(q_1_cur, q_2_cur));
    end
    
    Fd = F_fun(q_1_d, q_2_d);
    disp(['Fd: ' Fd]);
    Fd_fun = matlabFunction(Fd);

    d_star = fminunc(Fd_fun, 0);

    q_1_fun = matlabFunction(q_1_d);
    q_1_star = q_1_fun(d_star);

    q_2_fun = matlabFunction(q_2_d);
    q_2_star = q_2_fun(d_star);

    F_star = Fd_fun(d_star);
    
    disp(['q_1_star: ' num2str(q_1_star)]);
    disp(['q_2_star: ' num2str(q_2_star)]);
    disp(['F_star: ' num2str(F_star)]);
    scatter(q_1_star, q_2_star, 'r');
    plot([q_1_cur, q_1_star], [q_2_cur, q_2_star], '-k');
    text(q_1_star, q_2_star, ['CG Round ' num2str(i)]);
    
    q_1_pre = q_1_cur;
    q_2_pre = q_2_cur;
    q_1_cur = q_1_star;
    q_2_cur = q_2_star;
end
