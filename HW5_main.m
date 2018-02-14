function HW5_main
    clc
    x = [0 1 2 3];
    y = [0 1 1 1];
    %HW5_LSCF_parabolic_lsqcurvefit(x, y);
    HW5_LSCF_parabolic_fsolve(x, y);

return

function HW5_LSCF_parabolic_lsqcurvefit(feat, pred)
    x = feat;
    y = pred;
    hold on
    plot(x, y, 'ro');
    title('fitting using lsqcurvefit');
    
    f = @(c, xdata) c(1).*xdata.^2 + c(2).*xdata + c(3);
    c_0 = [1 1 1];
    
    [c,resnorm,~,exitflag,output] = lsqcurvefit(f,c_0,x,y);
    disp(c);
    plot(x, y, '.-', x, f(c,x), 'r:x');
    
return

function HW5_LSCF_parabolic_fsolve(feat, pred)
    x = feat;
    y = pred;
    hold on
    plot(x, y, 'ro');
    title('fitting using fminunc');

    %f = @(c, xdata) c(1).*xdata.^2 + c(2).*xdata + c(3);
    %Fsumsquares = @(c) sum((f(c, x) - y).^2);
    
    syms L M N
    f = L.*x.^2 + M.*x + N;
    Fsumsquares = sum((f(L, M, N, x) - y).^2);
    g = gradient(Fsumsquares, [L, M, N]);
    fn = {matlabFunction(g(1)), matlabFunction(g(2)), matlabFunction(g(3))};
    
    fn1 = @(a)[fn{1}(a(1), a(2), a(3)), fn{2}(a(1), a(2), a(3)), fn{3}(a(1), a(2), a(3))]; 
    
    %opts = optimoptions('fminunc','Algorithm','quasi-newton');
    a_0 = [1 1 1];
    
    a_min = fsolve(fn1, a_0);
    
    %[xunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,c_0,opts);
    %disp(xunc);
    %plot(x, y, '.-', x, f(xunc,x), 'r:x');
    
    disp(a_min);
    
return

function HW5_LSCF_exponential

return
