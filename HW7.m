clc

syms x y
b = (1 - x)^2 + (y - x^2)^2;
g = gradient(b, [x, y]);
fcontour(b, 'LevelStep', 2);
grid on;
axis([0 3 0 3]);
axis square
hold on

b_fun = matlabFunction(b);

% derivative
s1 = matlabFunction(g(1));
s2 = matlabFunction(g(2));

% initial value
x_cur = 2;
y_cur = 2;
b_cur = b_fun(y_cur, y_cur);

disp(['x_cur: ' num2str(x_cur)]);
disp(['y_cur: ' num2str(y_cur)]);
disp(['b_cur: ' num2str(b_cur)]);
scatter(x_cur, y_cur, 'b');
text(x_cur, y_cur, 'Initial point');

for i = 1:4
    syms d

    x_d = x_cur + d*(-1)*s1(x_cur, y_cur);
    y_d = y_cur + d*(-1)*s2(x_cur, y_cur);

    bd = (1 - x_d)^2 + (y_d - x_d^2)^2;
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
    text(x_star, y_star, ['Round ' num2str(i)]);
    
    x_cur = x_star;
    y_cur = y_star;
end

