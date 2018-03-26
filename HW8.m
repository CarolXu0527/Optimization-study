clc

syms x y
b = log(((x^2)*cos(y) + (y^2)*sin(x))^2 + (x^2)*(y^2) + (x+y)^2 + 1);
g = gradient(b, [x, y]);
fcontour(b, 'LevelStep', 2);
grid on;
axis([-1 5 -1 5]);
axis square
hold on

b_fun = matlabFunction(b);

% derivative
s1 = matlabFunction(g(1));
s2 = matlabFunction(g(2));

% initial value for SD
x_cur = 4;
y_cur = 1;
b_cur = b_fun(y_cur, y_cur);

disp(['x_cur: ' num2str(x_cur)]);
disp(['y_cur: ' num2str(y_cur)]);
disp(['b_cur: ' num2str(b_cur)]);
scatter(x_cur, y_cur, 'b');
text(x_cur, y_cur, 'Initial point');

for i = 1:6
    syms d

    x_d = x_cur + d*(-1)*s1(x_cur, y_cur);
    y_d = y_cur + d*(-1)*s2(x_cur, y_cur);

    bd = log(((x_d^2)*cos(y_d) + (y_d^2)*sin(x_d))^2 + (x_d^2)*(y_d^2) + (x_d+y_d)^2 + 1);
    bd_fun = matlabFunction(bd);

    d_star = fminunc(bd_fun, 0);

    x_fun = matlabFunction(x_d);
    x_star = x_fun(d_star);

    y_fun = matlabFunction(y_d);
    y_star = y_fun(d_star);

    b_star = bd_fun(d_star);
    
    disp(['x_star: ' num2str(x_star)]);
    disp(['y_star: ' num2str(y_star)]);
    disp(['b_star: ' num2str(b_star)]);
    scatter(x_star, y_star, 'b');
    plot([x_cur, x_star], [y_cur, y_star], '-k');
    text(x_star, y_star, ['SD Round ' num2str(i)]);
    
    x_cur = x_star;
    y_cur = y_star;
end

% initial value for CG
x_cur = 4;
y_cur = 1;
b_cur = b_fun(y_cur, y_cur);

disp(['x_cur: ' num2str(x_cur)]);
disp(['y_cur: ' num2str(y_cur)]);
disp(['b_cur: ' num2str(b_cur)]);
scatter(x_cur, y_cur, 'r');
text(x_cur, y_cur, 'Initial point');

for i = 1:6
    syms d

    if(i == 1)
        x_d = x_cur + d*(-1)*s1(x_cur, y_cur);
        y_d = y_cur + d*(-1)*s2(x_cur, y_cur);

    else
        tmp1 = s1(x_cur, y_cur)^2 + s2(x_cur, y_cur)^2;
        tmp2 = s1(x_pre, y_pre)^2 + s2(x_pre, y_pre)^2;
        x_d = x_cur + d*((-1)*s1(x_cur, y_cur) + (tmp1/tmp2)*s1(x_cur, y_cur));
        y_d = x_cur + d*((-1)*s2(x_cur, y_cur) + (tmp1/tmp2)*s2(x_cur, y_cur));
    end
    bd = log(((x_d^2)*cos(y_d) + (y_d^2)*sin(x_d))^2 + (x_d^2)*(y_d^2) + (x_d+y_d)^2 + 1);
    bd_fun = matlabFunction(bd);

    d_star = fminunc(bd_fun, 0);

    x_fun = matlabFunction(x_d);
    x_star = x_fun(d_star);

    y_fun = matlabFunction(y_d);
    y_star = y_fun(d_star);

    b_star = bd_fun(d_star);
    
    disp(['x_star: ' num2str(x_star)]);
    disp(['y_star: ' num2str(y_star)]);
    disp(['b_star: ' num2str(b_star)]);
    scatter(x_star, y_star, 'r');
    plot([x_cur, x_star], [y_cur, y_star], '-k');
    text(x_star, y_star, ['CG Round ' num2str(i)]);
    
    x_pre = x_cur;
    y_pre = y_cur;
    x_cur = x_star;
    y_cur = y_star;
end
