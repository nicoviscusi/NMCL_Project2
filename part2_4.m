clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

% Set to true if you want to see the animation of the solutions over time
animation = "False";

%% Resolution of the problem

% Definition of parameters
g = 1;
u = 0.25;

% Spatial domain
xspan = [0, 2];

% Temporal domain
tspan = [0, 2];

% Initial conditions
h0 = @(x) 1 - 0.1 * sin(pi * x);
m0 = @(x) 0 * x;

% Number of grid points
N = 200;

% Number of time steps
CFL = 0.5;

% Note that max(h0) = 1.5
% dt = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
% K = round((tspan(end) - tspan(1)) / dt);

% Source function
S = @(x, t) [0*x; 0*t];

% We integrate the source term exactly 
Sa = set_Sa(2);

% Here we use periodic boundary condition as the option ('peri')
bc = 'peri';

% Choose order for WENO reconstruction
k = 3;

% Solve the problem
[h, m, xc, tvec] = solver(xspan, tspan, N, ...
    CFL, h0, m0, @LaxFriedrichs, @flux_phys, Sa, bc, k);

%%
% We visualize the solution
if animation == "True"
    figure(1)
    for i = 1 : 20 : length(tvec)

        subplot(2, 1, 1)
        plot(xc, h(:, i), 'LineWidth', 2)
        hold on
   %     plot(xc, h0(xc - tvec(i)), '--', 'Linewidth', 2)
        %title(['$h(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            %'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.9 1.1]);
        hold off
        %legend('Numerical solution', 'Exact solution', ...
            %'Interpreter', 'latex')
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xc, m(:, i), 'LineWidth', 2)
        hold on
    %    plot(xc, u * h0(xc - tvec(i)), '--',  'Linewidth', 2)
        %title(['$m(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            %'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([-0.1 0.1])
        hold off
        %legend('Numerical solution', 'Exact solution', ...
            %'Interpreter', 'latex')
        set(gca, 'Fontsize', 20)
        drawnow

    end
end



%%  Error analysis 

% Generate a reference solution
[h1_ex, m1_ex, xvec1_ex, ~] = solver(xspan, tspan, 1000, ...
    CFL, h0, m0, @LaxFriedrichs, @flux_phys, Sa, bc, k);

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(3:6);

% Note that we cannot solve for smalle values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec = zeros(size(N_vec));
err_m_vec = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);
    % k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
    % K = round((tspan(end) - tspan(1)) / k);
    [h1, m1, xc, tvec] = solver(xspan, tspan, N, ...
        CFL, h0, m0, @LaxFriedrichs, @flux_phys, Sa, bc, k);

    % We now want to compare h1(:, end) with h1_ex(:, end), 
    % but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xvec1
    h1_interp = interp1(xc, h1(:, end), xvec1_ex);
    m1_interp = interp1(xc, m1(:, end), xvec1_ex);
    err_h_vec(i) =  1/sqrt(N) * norm(h1_interp' - h1_ex(:, end)); 
    err_m_vec(i) =  1/sqrt(N) * norm(m1_interp' - m1_ex(:, end));
end


% Plot the error
figure(2)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--", delta_x_vec, delta_x_vec.^3, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=0.5\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", "interpreter", ...
    "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on


subplot(2,1,2)
loglog(delta_x_vec, err_m_vec, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--", delta_x_vec, delta_x_vec.^3, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=0.5\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", "interpreter", ...
    "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)