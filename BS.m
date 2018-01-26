T = @(x) 1/7*sqrt(50^2 + x^2)+1/2*sqrt(50^2 + (100-x)^2);

X_0 = 0;
T_Pre = T(X_0);

X_Delta = 30;
X_Cur = X_0 + X_Delta;
T_Cur = T(X_Cur);

while T_Pre > T_Cur && X_Cur + X_Delta <= 100
    T_Pre = T_Cur;
    
    X_Cur = X_Cur + X_Delta;
    T_Cur = T(X_Cur);
end

fprintf('\n Current X before loop is = %.2f', X_Cur);

X_Delta_BS = 18;

while X_Delta_BS >= 1 || T(X_Cur+X_Delta_BS)/T_Cur > 0.01
    
        X_left = X_Cur - X_Delta_BS;
        T_tmp_1 = T(X_left);

        X_right = X_Cur + X_Delta_BS;
        T_tmp_2 = T(X_right);
    
    if(T_Cur < T_tmp_1 && T_Cur < T_tmp_2)
        %T_Cur = T_pre;
        fprintf('\n Type 1 Current X is = %.2f', X_Cur);
    end
    
    if(T_tmp_1 < T_tmp_2 && T_tmp_1 < T_Cur)
        T_Cur = T_tmp_1;
        X_Cur = X_left;
        fprintf('\n Type 2 Current X is = %.2f', X_left);
    end
    
    if(T_tmp_2 < T_tmp_1 && T_tmp_2 < T_Cur)
        T_Cur = T_tmp_2;
        X_Cur = X_right;
        fprintf('\n Type 3 Current X is = %.2f', X_right);
    end
    
    X_Delta_BS = X_Delta_BS/2;
end
 
fprintf('\n Final X is = %.2f', X_Cur);
fprintf('\n Final T is = %.2f', T(X_Cur));

fprintf('\n Min is = %.2f by minbnd', fminbnd(T,0,100));