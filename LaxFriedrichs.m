% Lax−Friedrichs flux at cell−interface

function fval = LaxFriedrichs(flux_phys, Ul, Ur, dx, dt)
% Separate, and compute component-wise

fval = 0.5 * (flux_phys(Ul) + flux_phys(Ur)) - 0.5 * dx/ dt * (Ur - Ul);
%fval = 0.5*(UL(2:end) + UR(1:end-1)) - 0.5*dx/dt*(UL(2:end)-UR(1:end-1));

return;






