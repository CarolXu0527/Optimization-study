function HW6_main
    clc;

    syms x y z
    F_syms = 2*x^2 + 2*y^2 + z^4 + 8*x + 4*y + 10*z + 2*exp(-(z/4)^2);
    fsurf(F_syms);
    grid on;
    axis([-4 4 -4 4]);
    axis square
    hold on

    x_0 = 4;
    y_0 = 8;
    z_0 = 4;
    delta = 1;
    range = [-1, 0 , 1];
    f = matlabFunction(F_syms);
    f_cur = f(x_0, y_0, z_0);

    f_min = Inf;
    f_max = -Inf;

    % 2 steps of marching grid search.
    % 1. find the direction
    % 2. find the minimum in the direction
    for i = 1:length(range)
        for j = 1:length(range)
            for k = 1:length(range)
                tmp = f(x_0 + delta*range(i), y_0 + delta*range(j), z_0 + delta*range(k));
                f_min = min(f_min, tmp);
                f_max = max(f_max, tmp);
                fprintf('\n Test on %.2f, %.2f, %.2f, %.2f', x_0 + delta*range(i), y_0 + delta*range(j), z_0 + delta*range(k), tmp);
                if (tmp < f_cur)
                    f_cur = tmp;
                    x_0_tmp = x_0 + delta*range(i);
                    y_0_tmp = y_0 + delta*range(j);
                    z_0_tmp = z_0 + delta*range(k);
                    fprintf('\n Found a lower point: %.2f, %.2f, %.2f, T is %.2f', x_0_tmp, y_0_tmp, z_0_tmp, f_cur);
                end
            end
        end
    end

    fprintf('\n *************************************************************************************');
    fprintf('\n Found a low point for initializing steepest descent: %.2f, %.2f, %.2f, T is %.2f', x_0_tmp, y_0_tmp, z_0_tmp, f_cur);

    % start to do steepest descent
    delta_x = Inf;
    delta_y = Inf;
    delta_z = Inf;
    delta_f = Inf;
    m = 0;

    g = gradient(F_syms, [x, y, z]);
    fn = {matlabFunction(g(1)), matlabFunction(g(2)), matlabFunction(g(3))};

    while((delta_x >= 0.01 || delta_y >= 0.01 || delta_z >= 0.01) && delta_f >= 0.01 && m < 15)
        s = -[fn{1}(x_0_tmp) fn{2}(y_0_tmp) fn{3}(z_0_tmp)];

        syms d
        xd = x_0_tmp + d*s(1);
        yd = y_0_tmp + d*s(2);
        zd = z_0_tmp + d*s(3);
        fd = f(xd, yd, zd);
        fd_fun = matlabFunction(fd);
        d_min = fminunc(fd_fun, [0]);
        x_star_fun = matlabFunction(xd);
        x_star = x_star_fun(d_min);
        y_star_fun = matlabFunction(yd);
        y_star = y_star_fun(d_min);
        z_star_fun = matlabFunction(zd);
        z_star = z_star_fun(d_min);
        f_star = fd_fun(d_min);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delta_x = abs(x_star - x_0_tmp);
        delta_y = abs(y_star - y_0_tmp);
        delta_z = abs(z_star - z_0_tmp);
        delta_f = abs(f_star - f_cur);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x_0_tmp = x_star;
        y_0_tmp = y_star;
        z_0_tmp = z_star;
        f_cur = f_star;
        m = m + 1;

        fprintf('\n %.2f, %.2f, %.2f, %.2f, %f', delta_x, delta_y, delta_z, delta_f, m);
        fprintf('\n Found a low point using steepest descent: %.2f, %.2f, %.2f, T is %.2f', x_star, y_star, z_star, f_star);
    end

end
