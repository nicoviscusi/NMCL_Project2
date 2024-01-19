function RHS = evalRHS(U, N, dt, dx, flux, flux_phys, bc, Crec, ...
    k, xf, dw, beta)

% EVALRHS - Evaluate the right-hand side (RHS) of a numerical scheme for a
%           hyperbolic conservation law.
%
%   RHS = EVALRHS(U, N, dt, dx, flux, flux_phys, bc, Crec, k, xf, dw, beta)
%   computes the RHS of a hyperbolic conservation law based on the given
%   solution vector 'U', numerical parameters, and flux functions.
%
%   Input:
%       U          - Solution matrix with two components.
%       N          - Number of grid points in the spatial domain.
%       dt         - Time step size.
%       dx         - Spatial grid spacing.
%       flux       - Flux function for the conservative variables.
%       flux_phys  - Physical flux.
%       bc         - String specifying the type of boundary condition.
%                    Supported values: 'peri' (periodic), 'open' 
%                    (open boundary).
%       Crec       - Reconstruction method coefficient 
%                    (for WENO reconstruction).
%       k          - Order of accuracy for WENO reconstruction.
%       xf         - Spatial grid locations.
%       dw         - WENO weights.
%       beta       - Parameter for WENO reconstruction.
%
%   Output:
%       RHS        - Computed right-hand side vector for 
%                    the conservation law.
%
%   Note:
%       This function applies appropriate boundary conditions,
%       performs WENO reconstruction, and computes the flux differences
%       to obtain the RHS.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% January 2024

% Apply appropriate boundary conditions
U = apply_bc(U, bc, k);

% Obtain reconstructed states
hl = zeros(1, N+2);
hr = zeros(1, N+2);

ml = zeros(1, N+2);
mr = zeros(1, N+2);

for i = 1:N+2
    [hl(i), hr(i)] = WENO(xf, U(1, i:(i+2*(k-1)))', k, Crec, dw, beta);
    [ml(i), mr(i)] = WENO(xf, U(2, i:(i+2*(k-1)))', k, Crec, dw, beta);
end

qr = [hr(2:N+1); mr(2:N+1)]; ql = [hl(2:N+1); ml(2:N+1)];
qm = [hr(1:N); mr(1:N)]; qp = [hl(3:N+2); ml(3:N+2)];

fluxval1 = flux(flux_phys, qr, qp, dx, dt);
fluxval2 = flux(flux_phys, qm, ql, dx, dt);

RHS = - (fluxval1-fluxval2);

end