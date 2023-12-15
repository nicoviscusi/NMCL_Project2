clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

% Set to true if you want to see the animation of the solutions over time
animation = "True";

%% Resolution of the problem

% Definition of parameters
g = 1;
u = 0.25;

% Spatial domain
xspan = [0, 2];

% Temporal domain
tspan = [0, 0.5];

% Initial conditions
h0 = @(x) 1 + 0.5 * sin(pi * x);
m0 = @(x) u * h0(x);

% Number of grid points
N = 1000;

% Number of time steps
CFL = 0.5;

% Note that max(h0) = 1.5
k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
K = round((tspan(end) - tspan(1)) / k);

% Source function
S = @(x, t) [pi/2 * (u - 1) * cos(pi * (x - t));
    pi/2 * cos(pi * (x - t)) * (- u + u^2 + g * h0(x - t))];

% Here we use periodic boundary condition as the option ('peri')
bc = 'peri';

% Solve the problem
[h, m, tvec, xvec, k, delta_x] = conservative_scheme(xspan, tspan, N, ...
    K, h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);


% We visualize the solution
if animation == "True"
    figure(1)
    for i = 1 : 20 : length(tvec)

        subplot(2, 1, 1)
        plot(xvec, h(:, i), 'LineWidth', 2)
        hold on
        plot(xvec, h0(xvec - tvec(i)), 'Linewidth', 2)
        title(['$h(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.4 1.6]);
        hold off
        legend('Numerical solution', 'Exact solution', ...
            'Interpreter', 'latex')
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec, m(:, i), 'LineWidth', 2)
        hold on
        plot(xvec, u * h0(xvec - tvec(i)), 'Linewidth', 2)
        title(['$m(x, t)$ at $t = $ ', num2str(tvec(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.1 0.4]);
        hold off
        legend('Numerical solution', 'Exact solution', ...
            'Interpreter', 'latex')
        set(gca, 'Fontsize', 20)
        drawnow

    end
end

%%  Error analysis 

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(4:10);

% Note that we cannot solve for smalle values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec = zeros(size(N_vec));
err_m_vec = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);
    k = CFL * (xspan(2) - xspan(1)) / N * 1 / (u + sqrt(g * 1.5));
    K = round((tspan(end) - tspan(1)) / k);
    T_f = 0.5;
    [h, m, ~, xvec, k, delta_x] = conservative_scheme(xspan, tspan, N, ...
        K, h0, m0, @lax_friedrichs_flux, @flux_phys, S, bc);
    err_h_vec(i) = 1/sqrt(N)*norm(h(:, end) -h0(xvec-T_f)');          
    err_m_vec(i) = 1/sqrt(N)*norm(m(:, end) - u*h0(xvec-T_f)');    
end


% Plot the error
figure(2)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=0.5\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on


subplot(2,1,2)
loglog(delta_x_vec, err_m_vec, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=0.5\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)