function RHS = evalRHS(U, N, dt, dx, flux, flux_phys, bc, Crec, ...
    k, xf, dw, beta)

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