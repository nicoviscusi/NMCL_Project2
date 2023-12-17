% Lax−Friedrichs flux at cell−interface

function fval = LaxFriedrichs(flux_phys, Ul, Ur, dx, dt)

fval = 0.5 * (flux_phys(Ul) + flux_phys(Ur)) - 0.5 * dx / dt * (Ur - Ul);

return;
