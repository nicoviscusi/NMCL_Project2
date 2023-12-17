function [h, m, xc, tvec] = solver(xspan, tspan, N, K, h0, m0, flux, flux_phys, Sa, bc, k)

dx      = (xspan(2) - xspan(1)) / N;
dt      = (tspan(2) - tspan(1)) / K;
tvec    = linspace(tspan(1), tspan(end), K + 1);
xc      = (xspan(1) + 0.5 * dx) : dx : (xspan(2) - 0.5 * dx);
xf      = linspace(xspan(1), xspan(2), N + 1);
h       = zeros(N, K + 1);
m       = zeros(N, K + 1);

% We compute the cell-averages of the initial condition for each variable
% We compute cell-averages of the initial condition
for j = 1 : N
    h(j, 1) = integral(h0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
    m(j, 1) = integral(m0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
end


% Initialize polynomial coefficients (of degree k-1)
Crec = zeros(k + 1, k);
for r=-1:k-1
    Crec(r+2,:) = ReconstructWeights(k,r);
end


% We can now start solving the problem

for i = 2 : length(tvec)

    hold = (h(:, i - 1))';
    hnew = hold;

    mold = (m(:, i - 1))';
    mnew = mold;

    % We use Runge-Kutta Strong Stability Preserving scheme to integrate in
    % time (SSP-RK3) - see exercise sheet for algorithm
    % Code mainly adapted from solutions

    % SSP-RK3 stage 1
    RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, Sa, bc, Crec, k, xf, tvec(i));
    hnew   = hold + dt * RHS(1);
    mnew   = mold + dt * RHS(2);

    % SSP-RK3 stage 2
    RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, Sa, bc, Crec, k, xf, tvec(i));
    hnew   = 3 * hold / 4.0 + (hnew + dt * RHS(1)) / 4.0;
    mnew   = 3 * mold / 4.0 + (mnew + dt * RHS(2)) / 4.0;

    % SSP-RK3 stage 3
    RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, Sa, bc, Crec, k, xf, tvec(i));
    hnew   = hold / 3.0 + 2.0 * (hnew + dt * RHS(1)) / 3.0;
    mnew   = mold / 3.0 + 2.0 * (mnew + dt * RHS(2)) / 3.0;

    h(:, i) = hnew';
    m(:, i) = mnew';

end
