% HW3 Several Methods for Single Variable Problems
clc

f = @(x) -x*exp(-x/5);
df = @(x) -exp(-x/5) + x*exp(-x/5)/5;


% Question 1
% Set the derivative equal to zero and solve for x* using ‘fzero’

ezplot(f,[0 10 -5 5]); % Plot function and derivative
hold on % Hold so marker can be added to show minimum
ezplot(df,[0 10 -5 5]); % Plot function and derivative
grid on

x_0 = 0;
x_min_zero = fzero(df, x_0);

plot(x_min_zero, f(x_min_zero), 'bo');

text(x_min_zero,fn(x_min_zero)-0.5,['x_{min} = ' num2str(x_min_zero)])

fprintf('x_min calculated by fzero is: %.2f\n', x_min_zero);

% Question 2
% Find the minimum using Matlab function ‘fminunc’

options = optimoptions('fminunc','GradObj','on');
fn = @obj_grad;

x_min_minunc = fminunc(fn, x_0, options);
fprintf('x_min calculated by fminunc is: %.2f\n', x_min_minunc);


% Find minimum using sequential quadratic approximation

% Use Monte Carlo / Quadratic hybrid method

% Use binary search / quadratic hybrid method