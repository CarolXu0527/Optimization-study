function HW5
    clc
    x = [0 1 2 3];
    y = [0 1 1 1];
    %HW5_LSCF_parabolic_lsqcurvefit(x, y);
    %HW5_LSCF_parabolic_fsolve(x, y);
    
    x_2 = [1 2 3 4];
    y_2 = [1 2 3 4];
    %HW5_LSCF_exponential_lsqcurvefit(x_2, y_2);
    HW5_LSCF_exponential_fsolve(x_2, y_2);

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
    title('fitting using fsolve');
    
    syms a b c
    
    Fsumsquares = 0;
    
    for i = 1:length(x)
        Fsumsquares = Fsumsquares + (y(i) - (a*x(i)^2 + b*x(i) + c))^2;
    end
    
    g = gradient(Fsumsquares, [a, b, c]);
    fn = {matlabFunction(g(1)), matlabFunction(g(2)), matlabFunction(g(3))};
    
    fn1 = @(a)[fn{1}(a(1), a(2), a(3)), fn{2}(a(1), a(2), a(3)), fn{3}(a(1), a(2), a(3))]; 
    
    %opts = optimoptions('fminunc','Algorithm','quasi-newton');
    a_0 = [1 1 1];
    
    a_min = fsolve(fn1, a_0);
    
    disp(a_min);
    plot(x, y, '.-', x, a_min(1)*x.^2 + x.*a_min(2) + a_min(3), 'r:x');
    
return

function HW5_LSCF_exponential_lsqcurvefit(feat, pred)
    x = feat;
    y = pred;
    hold on
    plot(x, y, 'ro');
    title('fitting using lsqcurvefit');
    
    f = @(c, xdata) c(1).*(exp(c(2).*xdata));
    c_0 = [1 1];
    
    c = lsqcurvefit(f,c_0,x,y);
    disp(c);
    plot(x, y, '.-', x, f(c,x), 'r:x');

return

function HW5_LSCF_exponential_fsolve(feat, pred)
    x = feat;
    y = pred;
    hold on
    plot(x, y, 'ro');
    title('fitting using fsolve');
    
    syms a b
    
    Fsumsquares = 0;
    
    for i = 1:length(x)
        Fsumsquares = Fsumsquares + (y(i) - (a.*(exp(b.*x(i)))))^2;
    end
    
    g = gradient(Fsumsquares, [a, b]);
    fn = {matlabFunction(g(1)), matlabFunction(g(2))};
    
    fn1 = @(a)[fn{1}(a(1), a(2)), fn{2}(a(1), a(2))]; 
    
    opts = optimoptions('fsolve','MaxFunctionEvaluations',400, 'MaxIterations', 400);
    a_0 = [1 1];
    
    a_min = fsolve(fn1, a_0, opts);
    
    disp(a_min);
    plot(x, y, '.-', x, a_min(1)*(exp(a_min(2).*x)), 'r:x');
    
return
