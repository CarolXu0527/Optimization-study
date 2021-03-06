clc
syms A B C

Vol = -A*B*C;
g = 2*B + 2*C - 2;
g_2 = A + 2*C - 1;
g_3 = -A;
g_4 = -B;
g_5 = -C;

% P = A + R*Hev(g)*g^2;
% R = 1;
%%% initiate the exterior penality function
P = Vol + heaviside(g)*(g^2) + heaviside(g_2)*(g_2^2) + heaviside(g_3)*(g_3^2) + heaviside(g_4)*(g_4^2) + heaviside(g_5)*(g_5^2);
P_g = gradient(P, [A, B, C]);

P_fun = matlabFunction(P);

% derivative
s_A = matlabFunction(P_g(1));
s_B = matlabFunction(P_g(2));
s_C = matlabFunction(P_g(3));

% initial value
A_cur = 0.3;
B_cur = 0.3;
C_cur = 0.1;
A_pre = Inf;
B_pre = Inf;
C_pre = Inf;

P_cur = P_fun(A_cur, B_cur, C_cur);
P_d_pre = Inf;
count = 1;

while ((abs(A_cur-A_pre) > 0.01) || (abs(B_cur-B_pre) > 0.01) || (abs(C_cur-C_pre) > 0.01)) && (abs(A_cur*B_cur*C_cur-A_pre*B_pre*C_pre) > 0.01) && (count < 26)
%%% convert variables from A & C to d
    disp('++++++++++++Start new iteration+++++++++++++++');
    syms d

    A_d = A_cur + d*(-1)*s_A(A_cur, B_cur, C_cur);
    A_fun = matlabFunction(A_d);
    B_d = B_cur + d*(-1)*s_B(A_cur, B_cur, C_cur);
    B_fun = matlabFunction(B_d);
    C_d = C_cur + d*(-1)*s_C(A_cur, B_cur, C_cur);
    C_fun = matlabFunction(C_d);

    P_d = -A_d*B_d*C_d + heaviside(2*B_d + 2*C_d - 2)*((2*B_d + 2*C_d - 2)^2) + heaviside(A_d + 2*C_d - 1)*((A_d + 2*C_d - 1)^2) + heaviside(-A_d)*((-A_d)^2) + heaviside(-B_d)*((-B_d)^2) + heaviside(-C_d)*((-C_d)^2);
    Pd_fun = matlabFunction(P_d);

    d_cur = 0;
    P_d_pre = Pd_fun(d_cur);
    d_delta = 1;
    d_cur = d_cur + d_delta;
    P_d_cur = Pd_fun(d_cur);

    %%% start to do binary search
    while P_d_pre > P_d_cur
        P_d_pre = P_d_cur;

        d_cur = d_cur + d_delta;
        P_d_cur = Pd_fun(d_cur);

    end

    d_cur = d_cur - d_delta;
    P_d_pre = Pd_fun(d_cur - d_delta);
    P_d_cur = Pd_fun(d_cur);

    disp(['P_d_pre: ' num2str(P_d_pre)]);
    disp(['P_d_cur: ' num2str(P_d_cur)]);
    
    while d_delta > 0.01
        d_delta = d_delta/2;
        P_d_pre = P_d_cur;
        
        d_left = d_cur - d_delta;
        Pd_tmp_1 = Pd_fun(d_left);

        d_right = d_cur + d_delta;
        Pd_tmp_2 = Pd_fun(d_right);

        if(P_d_cur < Pd_tmp_1 && P_d_cur < Pd_tmp_2)
            fprintf('Type 1 Current d is = %.2f\n', d_cur);
        end

        if(Pd_tmp_1 < Pd_tmp_2 && Pd_tmp_1 < P_d_cur)
            P_d_cur = Pd_tmp_1;
            d_cur = d_left;
            fprintf('Type 2 Current d is = %.2f\n', d_left);
        end

        if(Pd_tmp_2 < Pd_tmp_1 && Pd_tmp_2 < P_d_cur)
            P_d_cur = Pd_tmp_2;
            d_cur = d_right;
            fprintf('Type 3 Current d is = %.2f\n', d_right);
        end

        P_d_star = Pd_fun(d_cur);    
        P_d_cur = P_d_star;
        disp(['P_d_pre: ' num2str(P_d_pre)]);
        disp(['P_d_star: ' num2str(P_d_cur)]);
    end

    A_pre = A_cur;
    A_cur = A_fun(d_cur);
    B_pre = B_cur;
    B_cur = B_fun(d_cur);
    C_pre = C_cur;
    C_cur = C_fun(d_cur);
    count = count + 1;

    disp('+++++++++++++++++++++++++++');
    disp(['A_star: ' num2str(A_cur)]);
    disp(['B_star: ' num2str(B_cur)]);
    disp(['C_star: ' num2str(C_cur)]);
    disp(['Vol: ' num2str(A_cur*B_cur*C_cur)]);
    disp(['P_d_star: ' num2str(P_d_cur)]);
end
 
