function intS = integrate_source(xa, xb, t, PROBLEM)

u = 0.25;
g = 1;

if PROBLEM == 1

    int_x_S = @(xa,xb) [0.5*(u-1).*(sin(pi*(t-xa))-sin(pi*(t-xb)));
        1/8*sin(pi*(t-xb)).*(g*sin(pi*(t-xb)) -4*(g+(u-1)*u))- ...
        1/8*sin(pi*(t-xa)).*(g*sin(pi*(t-xa)) -4*(g+(u-1)*u))];
    intS = int_x_S(xa,xb);

elseif PROBLEM == 2

    int_x_S = @(xa,xb) [0*xa; 0*xa];
    intS = int_x_S(xa, xb);

end


return