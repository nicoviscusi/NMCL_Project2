function [h, m, xc, tvec] = solver(xspan, tspan, N, CFL, h0, m0, flux, flux_phys, Sa, bc, k)

g = 1;
dx      = (xspan(2) - xspan(1)) / N;
xc      = (xspan(1) + 0.5 * dx) : dx : (xspan(2) - 0.5 * dx);
xf      = linspace(xspan(1), xspan(2), N + 1);
h       = zeros(N,1);
m       = zeros(N,1);
tvec = [0];

% We compute the cell-averages of the initial condition for each variable
% We compute cell-averages of the initial condition
for j = 1 : N
    h(j) = integral(h0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
    m(j) = integral(m0, xf(j), xf(j+1), 'AbsTol', 1e-14) / dx;
end


% Initialize polynomial coefficients (of degree k-1)
% Crec = zeros(k + 1, k);
% for r=-1:k-1
%     Crec(r+2,:) = eval_crj(k,r);
% end
Crec = eval_crj(k);

% Store the initial condtion in q, to have size(q) = [2, Nspacenodes]
q = [h';m'];

% We can now start solving the problem
t = 0;
while (t<tspan(2))
   
    vel = q(2, :)./q(1,:);
    dt = dx*CFL/max(abs(vel) + sqrt(g*q(1,:)));

    % DEBUGGING
    if dt < 1e-5
        disp(h)
        disp(m)
        error("The solution is exploding")
    end
    
    if(t+dt>=tspan(2))
       dt = tspan(2) - t;
    end
    qold = q;
    Source = zeros(2,length(xf)-1);
    for j =1:length(xf)-1
        Source(:,j) = Sa(xf(j), xf(j+1), t, t+dt);
    end
    % SSP-RK3 Stage1
    RHS = evalRHS(q, N, dt, dx, flux, flux_phys, bc, Crec, k, xf, t);
    q = qold + (dt)*RHS + Source/dx;
    
    % SSP-RK3 Stage2
    RHS = evalRHS(q, N, dt, dx, flux, flux_phys, bc, Crec, k, xf, t);
    q = 3*qold/4.0 + (q + (dt)*RHS)/4.0+ Source/dx/4;

    % SSP-RK3 Stage3
    RHS = evalRHS(q, N, dt, dx, flux, flux_phys, bc, Crec, k, xf, t);
    q = qold/3.0 + 2.0*(q + (dt)*RHS)/3.0+ Source/dx*2/3;
    
    t = t + dt;
    
    % Store the new solution
    h = [h, q(1, :)'];
    m = [m, q(2, :)'];
    tvec = [tvec, t];

end




% for i = 2 : length(tvec)
% 
%     hold = (h(:, i - 1))';
%     hnew = hold;
% 
%     mold = (m(:, i - 1))';
%     mnew = mold;
% 
% 
%     % We use Runge-Kutta Strong Stability Preserving scheme to integrate in
%     % time (SSP-RK3) - see exercise sheet for algorithm
%     % Code mainly adapted from solutions
% 
%     Source = Sa(xf(1:end-1), xf(2:end), tvec(i), tvec(i) + dt) / (dx * dt);
% 
%     % SSP-RK3 stage 1
%     RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, bc, Crec, k, xf, tvec(i));
%     hnew   = hold + dt * (RHS(1,:) + Source(1, :));
%     mnew   = mold + dt * (RHS(2,:) + Source(2, :));
% 
%     Source = Sa(xf(1:end-1), xf(2:end), tvec(i), tvec(i) + 2*dt) / (dx * dt);
%     % SSP-RK3 stage 2
%     RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, bc, Crec, k, xf, tvec(i));
%     hnew   = 3 * hold / 4.0 + (hnew + dt * (RHS(1,:) + Source(1, :))) / 4.0 ;
%     mnew   = 3 * mold / 4.0 + (mnew + dt * (RHS(2,:) + Source(2, :))) / 4.0 ;
% 
%     Source = Sa(xf(1:end-1), xf(2:end), tvec(i), tvec(i) + dt+ dt/2) / (dx * dt);
%     % SSP-RK3 stage 3
%     RHS = evalRHS([hnew; mnew], N, dt, dx, flux, flux_phys, bc, Crec, k, xf, tvec(i));
%     hnew   = hold / 3.0 + 2.0 * (hnew + dt * (RHS(1,:) + Source(1, :))) / 3.0;
%     mnew   = mold / 3.0 + 2.0 * (mnew + dt * (RHS(2,:) + Source(2, :))) / 3.0;
% 
%     h(:, i) = hnew';
%     m(:, i) = mnew';
% 
% end
