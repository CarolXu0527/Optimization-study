% Minimize the surface area of a cylindrical can that holds at least 500 mL of liquid using both interior and exterior penalty function methods

% Use steepest descent with Matlab ‘fminunc’ function for 1-D search.  Initial design:  D=100mm, H=100mm

% Exit criteria is satisfied when change in f* < 1% or change in both design variables is < 1%.
clc
syms D H
A = pi*D*H + (pi/2)*D^2;
V = (D^2)*H*pi/4;
g = -(D^2)*H*pi/4 + 500;

% P = A + R*Hev(g)*g^2;
P = pi*D*H + (pi/2)*D^2 + heaviside(-(D^2)*H*pi/4 + 500)*(-(D^2)*H*pi/4 + 500)^2;
P_g = gradient(P, [D, H]);

fcontour(P);
grid on;

P_fun = matlabFunction(P);

% derivative
s_D = matlabFunction(P_g(1));
s_H = matlabFunction(P_g(2));

% initial value
D_cur = 100;
H_cur = 100;

P_cur = P_fun(D_cur, H_cur);






