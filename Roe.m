% Computes the Roe flux at cell−interfaces
function fval = Roe(flux_phys, Ul, Ur, dx, dt)

N    = length(Ul);
fval = zeros(2, N);

rUl  = sqrt(Ul(1,:));
rUr  = sqrt(Ur(1,:));

uavg = (Ul(2,:)./rUl+Ur(2,:)./rUr)./(rUl+rUr);
cavg = sqrt((Ul(1,:)+Ur(1,:))/2);

L1   = uavg-cavg;
L2   = uavg+cavg;
A11  = 0.5*(L2.*abs(L1)- L1.*abs(L2))./cavg;
A12  = 0.5*(abs(L2)- abs(L1))./cavg;
A21  = 0.5*L1.*L2.*(abs(L1)- abs(L2))./cavg;
A22  = 0.5*(L2.*abs(L2)- L1.*abs(L1))./cavg;
dU   = Ur-Ul;

fval = 0.5*(flux_phys(Ul) + flux_phys(Ur)) - 0.5*[A11.*dU(1,:) + A12.*dU(2,:); A21.*dU(1,:) + A22.*dU(2,:)];

return;

