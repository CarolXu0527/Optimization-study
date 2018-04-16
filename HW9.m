% Minimize the surface area of a cylindrical can that holds at least 500 mL of liquid using both interior and exterior penalty function methods

% Use steepest descent with Matlab ‘fminunc’ function for 1-D search.  Initial design:  D=100mm, H=100mm

% Exit criteria is satisfied when change in f* < 1% or change in both design variables is < 1%.
clc
syms D H
A = pi*D*H + (pi/2)*(D^2);
V = (D^2)*H*pi/4;
g = (D^2)*H*(-pi/4) + 500;
g_2 = -D;

% P = A + R*Hev(g)*g^2;
% R = 1;
P = A + heaviside(g)*(g^2) + heaviside(g_2)*(g_2^2);
P_g = gradient(P, [D, H]);

fcontour(P);
grid on;
axis([-10 100 -10 100]);
hold on

P_fun = matlabFunction(P);

% derivative
s_D = matlabFunction(P_g(1));
s_H = matlabFunction(P_g(2));

% initial value
D_cur = 100;
H_cur = 100;
D_pre = Inf;
H_pre = Inf;

P_cur = P_fun(D_cur, H_cur);
P_pre = Inf;

while ((abs((D_cur-D_pre)/D_cur) > 0.01) || (abs((H_cur-H_pre)/H_cur) > 0.01)) && (abs((P_cur-P_pre)/P_cur) > 0.01)
    syms d

    D_d = D_cur + d*(-1)*s_D(D_cur, H_cur);
    H_d = H_cur + d*(-1)*s_H(D_cur, H_cur);
    
    P_d = pi*D_d*H_d + (pi/2)*(D_d^2) + heaviside((D_d^2)*H_d*(-pi/4) + 500)*(((D_d^2)*H_d*(-pi/4) + 500)^2) + heaviside(-D_d)*((-D_d)^2);
    Pd_fun = matlabFunction(P_d);

    d_star = fminunc(Pd_fun, 0);
    
    D_fun = matlabFunction(D_d);
    D_star = D_fun(d_star);
    
    H_fun = matlabFunction(H_d);
    H_star = H_fun(d_star);
    
    P_star = Pd_fun(d_star);
    
    disp(['D_star: ' num2str(D_star)]);
    disp(['H_star: ' num2str(H_star)]);
    disp(['P_star: ' num2str(P_star)]);
    scatter(D_star, H_star, 'b');
    plot([D_cur, D_star], [H_cur, H_star], '-k');
    text(D_star, H_star, ['D = ' num2str(D_star) ' H = ' num2str(H_star)]);
    
    D_pre = D_cur;
    H_pre = H_cur;
    
    D_cur = D_star;
    H_cur = H_star;
    
    P_pre = P_cur;
    P_cur = P_star;
    
end
    
    
