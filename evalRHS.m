% Function evaluates the RHS

function RHS = evalRHS(U, N, dt, dx, flux, flux_phys, Sa, bc, Crec, k, xf, time)

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

% h_flux_val =  flux(flux_phys, hl(2:end), hr(1:end-1), dx, dt);
% m_flux_val =  flux(flux_phys, ml(2:end), mr(1:end-1), dx, dt);
% 
% RHS_h = -(h_flux_val(2:end) - h_flux_val(1:end-1));
% RHS_m = -(m_flux_val(2:end) - m_flux_val(1:end-1));
% 
% RHS = [RHS_h; RHS_m];
% RHS = RHS + Sa(xf(1:end-1),xf(2:end), time, time + dt) / (dx * dt);

%%%%%% Something has to be fixed here probably %%%%%%%%%%


% ql = [hl(2:N+1); ml(2:N+1)]; qr = [hr(2:N+1); mr(2:N+1)]; 
% qp = [hl(3:N+2); ml(3:N+2)]; qm = [hr(1:N); mr(1:N)];

UL = [hl; ml];
UR = [hr; mr];

% fluxval = flux(flux_phys, UR(:,1:end-1), UL(:,2:end), dx, dt);
fluxval = flux(flux_phys, UL(:,1:end-1), UR(:,2:end), dx, dt);
%fluxval = flux(flux_phys, UL, UR, dx, dt);

RHS = - 1 / dx * (fluxval(:,2:end) - fluxval(:,1:end-1)) ...
     + Sa(xf(1:end-1),xf(2:end), time, time + dt) / (dx*dt);

% RHS = - 1 / dx * (flux(flux_phys, qr, qp, dx, dt) - flux(flux_phys, qm, ql, dx, dt)) + ...
%     + Sa(xf(1:end-1), xf(2:end), time, time + dt)/dt/dx;

end