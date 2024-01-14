% Function evaluates the RHS

function RHS = evalRHS(U, N, dt, dx, flux, flux_phys, bc, Crec, k, xf, time)

% UL and UR are left and right reconstructed values for the cells
% i = 0,..., N+1
% Do not confuse them with the left and right values at interfaces

U = apply_bc(U, bc, k);

% Obtain reconstructed states
hl = zeros(1,N+2);
hr = zeros(1,N+2);

ml = zeros(1,N+2);
mr = zeros(1,N+2);

for i = 1:N+2
    [hl(i), hr(i)] = WENO(U(1, i:(i+2*(k-1)))', k, Crec);
    [ml(i), mr(i)] = WENO(U(2, i:(i+2*(k-1)))', k, Crec);
end

qr = [hr(2:N+1); mr(2:N+1)]; ql = [hl(2:N+1); ml(2:N+1)];
qm = [hr(1:N); mr(1:N)]; qp = [hl(3:N+2); ml(3:N+2)];

% Lax Friedrich flux - or something simple
% dq = - (MaxwellLF(qr,qp,ep,mu,k/h,maxvel) - ...
%              MaxwellLF(qm,ql,ep,mu,k/h,maxvel))/h;
fluxval1 = flux(flux_phys, qr, qp, dx, dt);
fluxval2 = flux(flux_phys, qm, ql, dx, dt);

RHS = - 1 / dx * (fluxval1-fluxval2);


% RHS = - 1 / dx * (flux(flux_phys, qr, qp, dx, dt) - flux(flux_phys, qm, ql, dx, dt)) + ...
%     + Sa(xf(1:end-1), xf(2:end), time, time + dt)/dt/dx;

end