function fval = LaxFriedrichs(flux_phys, Ul, Ur, dx, dt)
% LAXFRIEDRICHS - Compute the Lax-Friedrichs flux at a cell interface.
%
%   fval = LAXFRIEDRICHS(flux_phys, Ul, Ur, dx, dt) computes the 
%          Lax-Friedrichs flux at a cell interface based on the
%          physical flux function, left-state 'Ul', right-state 'Ur',
%          spatial grid spacing 'dx', and time step size 'dt'.
%
%   Input:
%       flux_phys - Physical flux function.
%       Ul        - Left-state values.
%       Ur        - Right-state values.
%       dx        - Spatial grid spacing.
%       dt        - Time step size.
%
%   Output:
%       fval      - Computed Lax-Friedrichs flux at the cell interface.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% January 2024


fval = 0.5 * (flux_phys(Ul) + flux_phys(Ur)) - 0.5 * dx/ dt * (Ur - Ul);

return;






