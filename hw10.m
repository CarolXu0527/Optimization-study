clc
syms A C

Vol = -A*A*C;
g = A + C - 1;
g_2 = 2*A + C -1;

% P = A + R*Hev(g)*g^2;
% R = 1;
%%% initiate the exterior penality function
P = A + heaviside(g)*(g^2) + heaviside(g_2)*(g_2^2);
P_g = gradient(P, [A, C]);

fcontour(P);
grid on;
axis([-5 20 -5 20]);
hold on

P_fun = matlabFunction(P);

% derivative
s_A = matlabFunction(P_g(1));
s_C = matlabFunction(P_g(2));

% initial value
A_cur = 0.1;
C_cur = 0.1;
A_pre = Inf;
C_pre = Inf;

P_cur = P_fun(A_cur, C_cur);
P_d_pre = Inf;

%%% convert variables from A & C to d
syms d

A_d = A_cur + d*(-1)*s_A(A_cur, C_cur);
C_d = C_cur + d*(-1)*s_C(A_cur, C_cur);
    
P_d = -A_d*A_d*C_d + heaviside(A_d + C_d - 1)*((A_d + C_d - 1)^2) + heaviside(2*A_d + C_d - 1)*((2*A_d + C_d - 1)^2);
Pd_fun = matlabFunction(P_d);

d_cur = 0;
P_d_pre = Pd_fun(d_cur);
fprintf('\n P_d_cur is = %.8f', P_d_cur);

d_delta = 0.1;
d_cur = d_cur + d_delta;
P_d_cur = Pd_fun(d_cur);

%%% start to do binary search
while P_d_pre > P_d_cur
    P_d_pre = P_d_cur;
    
    d_cur = d_cur + d_delta;
    P_d_cur = Pd_fun(d_cur);
    
end

count = 1;

while ((abs((A_cur-A_pre)/A_cur) > 0.01) || (abs((C_cur-C_pre)/C_cur) > 0.01)) && (abs((P_d_cur-P_d_pre)/P_d_cur) > 0.01) && (count < 26)
    d_left = d_cur - d_delta;
    Pd_tmp_1 = Pd_fun(d_left);

    d_right = d_cur + d_delta;
    Pd_tmp_2 = Pd_fun(d_right);
    
    count = count + 1;
    if(P_d_cur < Pd_tmp_1 && P_d_cur < Pd_tmp_2)
        %T_Cur = T_pre;
        fprintf('\n Type 1 Current X is = %.2f', d_cur);
    end
    
    if(Pd_tmp_1 < Pd_tmp_2 && Pd_tmp_1 < P_d_cur)
        P_d_cur = Pd_tmp_1;
        d_cur = d_left;
        fprintf('\n Type 2 Current X is = %.2f', d_left);
    end
    
    if(Pd_tmp_2 < Pd_tmp_1 && Pd_tmp_2 < P_d_cur)
        P_d_cur = Pd_tmp_2;
        d_cur = d_right;
        fprintf('\n Type 3 Current X is = %.2f', d_right);
    end
    
    d_delta = d_delta/2;
    
    A_fun = matlabFunction(A_d);
    A_star = A_fun(d_cur);
    
    C_fun = matlabFunction(C_d);
    C_star = C_fun(d_cur);
    
    P_d_star = Pd_fun(d_cur);
    
    A_pre = A_cur;
    C_pre = C_cur;
    
    A_cur = A_star;
    C_cur = C_star;
    
    P_d_pre = P_d_cur;
    P_d_cur = P_d_star;
end

disp(['A_star: ' num2str(A_cur)]);
disp(['C_star: ' num2str(C_cur)]);
disp(['P_d_star: ' num2str(P_d_cur)]);










    
