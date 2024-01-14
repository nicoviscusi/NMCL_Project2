function Sa = set_Sa(PROBLEM)
g = 1;

if PROBLEM == 1
    % Using exact integrated expression for the source terms (integrated in
    % space and time.
    intCos_pixmt =@(xa,xb,ta,tb) (sin(pi*xb) - sin(pi*xa)).*(sin(pi*tb) - sin(pi*ta))/pi^2 ...
        +(cos(pi*xb) - cos(pi*xa)).*(cos(pi*tb) - cos(pi*ta))/pi^2;
    intSin_2pixmt =@(xa,xb,ta,tb) (cos(2*pi*xb) - cos(2*pi*xa)).*(sin(2*pi*ta) ...
        - sin(2*pi*tb))/(2*pi)^2 ...
        -(sin(2*pi*xb) - sin(2*pi*xa)).*(cos(2*pi*ta) ...
        - cos(2*pi*tb))/(2*pi)^2;

    Sa = @(xa,xb,ta,tb) [0.5*pi*(0.25-1)*intCos_pixmt(xa,xb,ta,tb) ;
        0.5*pi*(0.25^2 -0.25+g)*intCos_pixmt(xa,xb,ta,tb) ...
        + g*pi*intSin_2pixmt(xa,xb,ta,tb)/8.0];
else
    Sa = @(xa,xb,ta,tb) [0*xa;
        0*xa];
end

return

