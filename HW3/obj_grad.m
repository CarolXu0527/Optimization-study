function [f, df] = obj_grad(x)
    f = -x*exp(-x/5);
    df = -exp(-x/5) + x*exp(-x/5)/5;
return