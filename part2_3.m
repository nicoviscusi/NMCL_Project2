clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

% Set to true if you want to see the animation of the solutions over time
animation = "True";

%% Resolution of the problem

% First problem (for integration of source term in space)
PROBLEM = 1;

% Definition of parameters
g = 1;
u = 0.25;

% Spatial domain
xspan = [0, 2];

% Temporal domain
tspan = [0, 1];

% Initial conditions
h0 = @(x) 1 + 0.5 * sin(pi * x);
m0 = @(x) u * h0(x);

% Number of grid points
N = 500;

% CFL condition
CFL = 0.5;

% Here we use periodic boundary condition as the option ('peri')
bc = 'peri';

% Choose order for WENO reconstruction
k = 2;

% Solve the problem
[h, m, xc, tvec] = solver(xspan, tspan, N, ...
    CFL, g, h0, m0, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

% We visualize the solution
if animation == "True"
    figure(1)
    for i = 1 : 1 : length(tvec)

        subplot(2, 1, 1)
        plot(xc, h(:, i), 'LineWidth', 2)
        hold on
        plot(xc, h0(xc - tvec(i)), '--', 'Linewidth', 2)
        title(['$h(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.5 1.5]);
        hold off
        legend('Numerical solution', 'Exact solution', ...
            'Interpreter', 'latex', 'Location', 'best')
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xc, m(:, i), 'LineWidth', 2)
        hold on
        plot(xc, u * h0(xc - tvec(i)), '--',  'Linewidth', 2)
        title(['$m(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0 0.4]);
        hold off
        legend('Numerical solution', 'Exact solution', ...
            'Interpreter', 'latex', 'Location', 'best')
        set(gca, 'Fontsize', 20)
        drawnow

    end
end



%%  Error analysis 

% We solve the same problem for different values of \Delta x
delta_x_vec = 2.^-(6:9);

N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec = zeros(size(N_vec));
err_m_vec = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);

    [h, m, xc, tvec] = solver(xspan, tspan, N, ...
        CFL, g, h0, m0, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);
    
    err_h_vec(i) = 1/sqrt(N) * norm(h(:, end) - h0(xc-tspan(2))');          
    err_m_vec(i) = 1/sqrt(N) * norm(m(:, end) - u*h0(xc-tspan(2))');    
end


% Plot the error in loglog 
figure(2)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^2, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^3, "--", "Linewidth", 2)
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=1\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", ...
    "interpreter", "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on


subplot(2,1,2)
loglog(delta_x_vec, err_m_vec, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^2, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^3, "--", "Linewidth", 2)
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=1\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", ...
    "interpreter", "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)