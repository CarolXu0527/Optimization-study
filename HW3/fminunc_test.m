% HW3 Several Methods for Single Variable Problems
% Question 2
% Find the minimum using Matlab function ‘fminunc'


clc
f = @(x) -x*exp(-x/5);
df = @(x) -exp(-x/5) + x*exp(-x/5)/5;
x_0 = 0;
options = optimoptions('fminunc','GradObj','on');
fn = [f, df];
x_min_minunc = fminunc(f, x_0);

fprintf('x_min calculated by fminunc is: %.2f\n', x_min_minunc);
