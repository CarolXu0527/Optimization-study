T = @(x, y) 1/7*sqrt(50^2 + x^2) + 1/4*sqrt(20^2 + (y-x)^2) + 1/2*sqrt(30^2 + (100-y)^2);

%%%%%%%%%%% Initialize marching grid %%%%%%%%%%%%%%
x_1 = 50;
x_2 = 50;
delta = 30;
range = [-1, 0, 1];
T_cur = T(x_1, x_2);

fprintf('\n Print initial T: %.2f', T_cur);
fprintf('\n Print initial Delta: %.2f', delta);

%%%%%%%%%%%%  ?x < 0.5m Change in T < 0.005 %%%%%%%%
while(delta >= 0.5 || T_delta >= 0.005)
    T_pre = Inf;
    T_delta = Inf;

    while(T_pre ~= T_cur)
        fprintf('\n ***************start new point*********************\n');
        fprintf('\n Current point: %.2f, %.2f, T is %.2f', x_1, x_2, T_cur)
        T_min = Inf;
        T_max = -Inf;
        T_pre = T_cur;

        for i = 1:length(range)
            for j = 1:length(range)
                tmp = T(x_1 + delta*range(i), x_2 + delta*range(j));
                T_min = min(T_min, tmp);
                T_max = max(T_max, tmp);
                fprintf('\n Test on %.2f, %.2f, %.2f', x_1 + delta*range(i), x_2 + delta*range(j), tmp);
                if (tmp < T_cur)   
                    T_cur = tmp;
                    x_1_tmp = x_1 + delta*range(i);
                    x_2_tmp = x_2 + delta*range(j);
                    fprintf('\n Found a lower point: %.2f, %.2f, T is %.2f', x_1_tmp, x_2_tmp, T_cur);
                end
            end
        end
        
        x_1 = x_1_tmp;
        x_2 = x_2_tmp;
        T_delta = T_max - T_min;
        fprintf('\n Current T_Delta: %.2f', T_delta);
        fprintf('\n Lowest point of current grid: %.2f, %.2f', x_1, x_2);    
    end
    
    %%%%%% Location of minimumm hasn't changed, so the grid size is
    %%%%%% reduced. 
    delta = delta/2;
    fprintf('\n Reduce delta to: %.2f', delta);
    fprintf('\n **************all**********************\n');
end

fprintf('\n Jump out with Delta: %.4f, T_delta %.4f', delta, T_delta);
fprintf('\n ***************************************\n');
fprintf('\n Pring final result: %.2f: ', T_cur);
