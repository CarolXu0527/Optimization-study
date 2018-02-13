
function HW3_main
clc


[x_ans_1, y_ans_1] = HW3_1_fzero;
disp('------- fzero approach ----------');
disp(['x = ' num2str(x_ans_1)]);
disp(['y = ' num2str(y_ans_1)]);

[x_ans_2, y_ans_2] = HW3_2_fminunc;
disp('------- fminunc approach ---------');
disp(['x = ' num2str(x_ans_2)]);
disp(['y = ' num2str(y_ans_2)]);

[x_ans_3, y_ans_3] = HW3_3_approximation(0);
disp('---- approximation approach ------');
disp(['x = ' num2str(x_ans_3)]);
disp(['y = ' num2str(y_ans_3)]);

[x_ans_4, y_ans_4] = HW3_4_hybrid_Monte_Carlo;
disp('----- hybrid with Monte Carlo ------');
disp(['x = ' num2str(x_ans_4)]);
disp(['y = ' num2str(y_ans_4)]);

[x_ans_5, y_ans_5] = HW3_5_hybrid_binary_search;
disp('----- hybrid with binary search ----');
disp(['x = ' num2str(x_ans_5)]);
disp(['y = ' num2str(y_ans_5)]);
return

% Set the derivative equal to zero and solve for x* using ‘fzero’
function [min_x, min_y] = HW3_1_fzero

f = @(x) -x.*exp(-x./5);
df = @(x) -exp(-x./5) + x.*exp(-x./5)/5;

ezplot(f, [-15 15 -15 15]);
hold on
grid on
axis square

x_0 = 0;
x_root = fzero(df, x_0);
plot(x_root, f(x_root), 'bo');

text(x_root + 0.5, f(x_root), ['x = ' num2str(x_root) ' y = ' num2str(f(x_root))])

min_x = x_root;
min_y = f(x_root);

return

% Find the minimum using Matlab function ‘fminunc’
function [min_x, min_y] = HW3_2_fminunc

x_0 = 0;
options = optimoptions('fminunc','GradObj','on','Algorithm','quasi-newton');

fn = @obj_grad;
x_min_minunc = fminunc(fn, x_0, options);

min_x = x_min_minunc;
min_y = obj_grad(x_min_minunc);

return

% Object used in solution 2
function [f, df] = obj_grad(x)
    f = -x*exp(-x/5);
    df = -exp(-x/5) + x*exp(-x/5)/5;
return

% Find minimum using sequential quadratic approximation
function [min_x, min_y] = HW3_3_approximation(x)

disp('****************************************************');
disp(['Starting point of SQA is: ' x]);
x_0 = x;
x_Delta = 1;
x_Delta_star_min = 0.05;
f_Delta_star_min = 0.05;
x_Delta_star = realmax;
f_Delta_star = realmax;

f = @(x) -x.*exp(-x./5);

x_1 = x_0 + x_Delta;
x_2 = x_0 + 2*x_Delta;

f_0 = f(x_0);
f_1 = f(x_1);
f_2 = f(x_2);

a = [x_0^2, x_0, 1; x_1^2, x_1, 1; x_2^2, x_2, 1];
b = [f_0; f_1; f_2];

c = a\b;

p = @(x) c(1)*x.^2 + x.*c(2) + c(3);
figure 
subplot(1, 2, 1)
fplot(f, [0 30], '-');
title('Original Function');
axis square

subplot(1, 2, 2)
fplot(p, [0 10], '-.*');
hold on
grid on
title('Sequential Quadratic Approximation');
axis square

x_star = -c(2)/(2*c(1));

while c(1) > 0 && (x_Delta_star >= x_Delta_star_min || f_Delta_star >= f_Delta_star_min)
    
    
    x_0 = x_star - x_Delta;
    x_1 = x_star;
    x_2 = x_star + x_Delta;
    
    f_0 = f(x_0);
    f_1 = f(x_1);
    f_2 = f(x_2);
    
    a = [x_0^2, x_0, 1; x_1^2, x_1, 1; x_2^2, x_2, 1];
    b = [f_0; f_1; f_2];
    
    c = a\b;
    
    p = @(x) c(1)*(x.^2) + x.*c(2) + c(3);
    x_Delta_star = abs(x_star - (-c(2)/(2*c(1))));
    f_Delta_star = abs(f(x_star)- f(-c(2)/(2*c(1))));
    x_star = -c(2)/(2*c(1));
    fplot(p, [0 10],'-.*');
    
    disp(['c(1): ' num2str(c(1))]);
    disp(['x_Delta: ' num2str(x_Delta)]);
    disp(['x_Delta_star: ' num2str(x_Delta_star)]);
    disp(['f_Delta_star: ' num2str(f_Delta_star)]);
    disp(['x_star: ' num2str(x_star)]);
    disp('Condition status after this round is shown as below:');
    disp(c(1) > 0);
    disp(x_Delta_star >= x_Delta_star_min);
    disp(f(x_star) >= f_Delta_star_min);
    
    if x_Delta_star >= x_Delta
        x_Delta = x_Delta/2;
        disp('hit delta_x changing point');
    end
    
    disp('****************************************************');
end

hold off

min_x = x_star;
min_y = f(x_star);

return

% Use Monte Carlo / Quadratic hybrid method
function [min_x, min_y] = HW3_4_hybrid_Monte_Carlo

d = 10;
x = random('uniform', 0, d, 5, 1);
f = -x.*exp(-x./5);

fplot(@(x) -x.*exp(-x./5), [0 10]);

hold on
grid on
axis square

[F, I] = sort(f);
disp(['x_star of MC: ' num2str(x(I(1))) ])
disp(['f_star of MC: ' num2str(F(1)) ])
plot(x(I(1)), F(1), 'rs');

[min_x, min_y] = HW3_3_approximation(x(I(1)));

return

% Use binary search / quadratic hybrid method
function [min_x, min_y] = HW3_5_hybrid_binary_search

f = @(x) -x.*exp(-x./5);

x_Cur = 0;
f_Pre = realmax;
f_Cur = f(x_Cur);
x_Delta = 3;

while f_Pre > f_Cur && x_Cur + x_Delta <= 30
    f_Pre = f_Cur;
    x_Cur = x_Cur + x_Delta;
    f_Cur = f(x_Cur);
end

x_0 = x_Cur - x_Delta;
[min_x, min_y] = HW3_3_approximation(x_0);

return
