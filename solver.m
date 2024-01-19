function [h, m, xc, tvec] = solver(xspan, tspan, N, CFL, g, h0, m0, ...
    flux, flux_phys, bc, k, PROBLEM)
% SOLVER - Solve a hyperbolic conservation law using the SSP-RK3 scheme.
%
%   [h, m, xc, tvec] = SOLVER(xspan, tspan, N, CFL, g, h0, m0, ...
%                             flux, flux_phys, bc, k, PROBLEM)
%   solves a hyperbolic conservation law using the SSP-RK3 scheme 
%   and returns the solutions for water height 'h', discharge 'm', 
%   spatial grid locations 'xc' and the time vector 'tvec'.
%
%   Input:
%       xspan     - Spatial domain [x_start, x_end].
%       tspan     - Time domain [t_start, t_end].
%       N         - Number of spatial grid points.
%       CFL       - Courant-Friedrichs-Lewy (CFL) number for adaptive 
%                   time-stepping.
%       g         - Gravitational constant (g=1).
%       h0        - Function handle for the initial water height
%                   condition.
%       m0        - Function handle for the initial discharge condition.
%       flux      - Flux function for the conservative variables.
%       flux_phys - Flux function in physical variables.
%       bc        - String specifying the type of boundary condition.
%                   Supported values: 'peri' (periodic), 'open'
%                   (open boundary).
%       k         - Order of accuracy for WENO reconstruction.
%       PROBLEM   - Integer specifying the problem type.
%                   Supported values: 1 or 2.
%
%   Output:
%       h         - Matrix of water height solutions over time.
%       m         - Matrix of discharge solutions over time.
%       xc        - Spatial grid locations.
%       tvec      - Time vector.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% January 2024


% Define preliminary variables
dx      = (xspan(2) - xspan(1)) / N;
xc      = (xspan(1) + 0.5 * dx) : dx : (xspan(2) - 0.5 * dx);
xf      = linspace(xspan(1), xspan(2), N + 1);
h       = zeros(N,1);
m       = zeros(N,1);
tvec    = 0;

% We compute cell-averages of the initial condition
for j = 1 : N
    h(j) = integral(h0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
    m(j) = integral(m0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
end

% Initialize polynomial coefficients (of degree k-1)
Crec = zeros(k + 1, k);
for r=-1:k-1
    Crec(r+2,:) = ReconstructWeights(k,r);
end

% Initialize linear weights
dw = LinearWeights(k,0);

% Compute smoothness indicator matrices
beta = zeros(k,k,k);
for r=0:k-1
    xl = -1/2 + [-r:1:k-r];
    beta(:,:,r+1) = betarcalc(xl,k);
end

% Store the initial condtion in q, to have size(q) = [2, Nspacenodes]
q = [h';m'];

% We can now start solving the problem
t = 0;

while (t < tspan(2))
   
    vel = q(2, :) ./ q(1,:);

    % Adaptive time-step
    dt = dx * CFL / max(abs(vel) + sqrt(g * q(1,:)));

    % DEBUGGING
    if dt < 1e-5
        disp(h)
        disp(m)
        error("The solution is exploding")
    end
    
    if(t + dt >= tspan(2))
       dt = tspan(2) - t;
    end

    qold = q;
    
    % Define variables that will store the source terms
    Source1 = zeros(2, length(xf)-1);
    Source2 = zeros(2, length(xf)-1);
    Source3 = zeros(2, length(xf)-1);

    % Integrate once for all the source term and evaluate it at different
    % time steps for later use from SSP-RK3
    if PROBLEM == 1

        Source1 = integrate_source(xf(1:end-1), xf(2:end), ...
            t, PROBLEM);
        Source2 = integrate_source(xf(1:end-1), xf(2:end), ...
            t + dt, PROBLEM);
        Source3 = integrate_source(xf(1:end-1), xf(2:end), ...
            t + 0.5 * dt, PROBLEM);

    end

 
    % We use Runge-Kutta Strong Stability Preserving scheme to integrate 
    % in time (SSP-RK3) - see exercise sheet for algorithm

    % SSP-RK3 Stage1
    RHS = evalRHS(q, N, dt, dx, flux, flux_phys, bc, Crec, ...
        k, xf, dw, beta);
    q1 = qold + dt/dx * (RHS + Source1);
    
    % SSP-RK3 Stage2
    RHS = evalRHS(q1, N, dt, dx, flux, flux_phys, bc, Crec, ...
        k, xf, dw, beta);
    q2 = 3*qold/4.0 + (q1 + dt/dx * (RHS + Source2))/4.0;

    % SSP-RK3 Stage3
    RHS = evalRHS(q2, N, dt, dx, flux, flux_phys, bc, Crec, ...
        k, xf, dw, beta);
    q3 = qold/3.0 + 2.0*(q2 + dt/dx * (RHS + Source3))/3.0;
    
    q = q3;

    t = t + dt;
    
    % Store the new solution
    h = [h, q(1, :)'];
    m = [m, q(2, :)'];
    tvec = [tvec, t];

end