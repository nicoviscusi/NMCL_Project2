function intS = integrate_source(xa, xb, t, PROBLEM)

% INTEGRATE_SOURCE - Compute the spatial integral of the source term for
%                    a given hyperbolic problem at a time t, using the 
%                    exact solution of the integral.
%
%   intS = INTEGRATE_SOURCE(xa, xb, t, PROBLEM) computes the spatial
%   integral of the source term for a specified problem at a given time
%   't'.
%
%   Input:
%       xa      - Left spatial boundary.
%       xb      - Right spatial boundary.
%       t       - Time at which the source term is evaluated.
%       PROBLEM - Integer specifying the problem type.
%                 Supported values: 1 or 2.
%
%   Output:
%       intS    - Computed spatial integral of the source term.
%
%   Parameters:
%       u       - Constant parameter (problem-dependent).
%       g       - Constant parameter (problem-dependent).
%
%   Problem Descriptions:
%       PROBLEM = 1: Computes the integral for the source term.
%       PROBLEM = 2: Returns zero for the source term.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% January 2024

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