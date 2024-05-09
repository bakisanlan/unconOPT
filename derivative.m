function dx = derivative(x,f)

eps = 10e-3;
x1 = x;
x2 = x+eps;


dx = ( f(x2) - f(x1) ) / (x2 - x1);

end