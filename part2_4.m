clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

% Set to true if you want to see the animation of the solutions over time
animation = "True";

%% Resolution of the problem (first set of initial conditions)

% Second problem
PROBLEM = 2;

% Definition of parameters
g = 1;
u = 0.25;

% Spatial domain
xspan = [0, 2];

% Temporal domain
tspan = [0, 1];

% Initial conditions
h01 = @(x) 1 - 0.1 * sin(pi * x);
m01 = @(x) 0 * x;

% Number of grid points
N = 500;

% Number of time steps
CFL = 0.5;

% Here we use periodic boundary condition as the option ('peri')
bc = 'peri';

% Choose order for WENO reconstruction
k = 2;

% Solve the problem
[h1, m1, xc1, tvec1] = solver(xspan, tspan, N, ...
    CFL, g, h01, m01, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

% We visualize the solution
if animation == "True"

    figure(1)

    for i = 1 : 1 : length(tvec1)

        subplot(2, 1, 1)
        plot(xc1, h1(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $ ', num2str(tvec1(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.9 1.1]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xc1, m1(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $ ', num2str(tvec1(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([-0.1 0.1])
        set(gca, 'Fontsize', 20)
        drawnow

    end

end



%%  Error analysis (first set of initial conditions)

% Generate a reference solution using sufficiently fine mesh
[h1_ex, m1_ex, xvec1_ex, ~] = solver(xspan, tspan, 1500, ...
    CFL, g, h01, m01, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(6:9);

N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec;
err_h_vec1 = zeros(size(N_vec));
err_m_vec1 = zeros(size(N_vec));

for i = 1 : length(N_vec)

    N = N_vec(i);

    [h1, m1, xc1, ~] = solver(xspan, tspan, N, ...
        CFL, g, h01, m01, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

    % We now want to compare h1(:, end) with h1_ex(:, end),
    % but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xc1
    h1ex_interp = interp1(xvec1_ex, h1_ex(:, end), xc1);
    m1ex_interp = interp1(xvec1_ex, m1_ex(:, end), xc1);

    % Compute norm 2 of the error
    err_h_vec1(i) =  1/sqrt(N) * norm(h1(:, end)' - h1ex_interp);
    err_m_vec1(i) =  1/sqrt(N) * norm(m1(:, end)' - m1ex_interp);

end


% Plot the error in loglog
figure(2)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec1, "o-", "Linewidth", 2)
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
loglog(delta_x_vec, err_m_vec1, "o-", "Linewidth", 2)
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



%% Resolution of the problem (second set of initial conditions)

% Initial conditions
h02 = @(x) 1 - 0.2 * sin(2 * pi * x);
m02 = @(x) 0.5 + 0*x;

% Number of grid points
N = 500;

% Solve the problem
[h2, m2, xc2, tvec2] = solver(xspan, tspan, N, ...
    CFL, g, h02, m02, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

% We visualize the solution
if animation == "True"

    figure(3)

    for i = 1 : 1 : length(tvec2)

        subplot(2, 1, 1)
        plot(xc2, h2(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $ ', num2str(tvec2(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.8 1.2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xc2, m2(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $ ', num2str(tvec2(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$m(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.4 0.6])
        set(gca, 'Fontsize', 20)
        drawnow

    end

end



%%  Error analysis (second set of initial conditions)

% Generate a reference solution
[h2_ex, m2_ex, xvec2_ex, ~] = solver(xspan, tspan, 1500, ...
    CFL, g, h02, m02, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(6:9);

N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec;
err_h_vec2 = zeros(size(N_vec));
err_m_vec2 = zeros(size(N_vec));

for i = 1 : length(N_vec)

    N = N_vec(i);

    [h2, m2, xc2, ~] = solver(xspan, tspan, N, ...
        CFL, g, h02, m02, @LaxFriedrichs, @flux_phys, bc, k, PROBLEM);

    % We now want to compare h1(:, end) with h1_ex(:, end),
    % but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xc2
    h2ex_interp = interp1(xvec2_ex, h2_ex(:, end), xc2);
    m2ex_interp = interp1(xvec2_ex, m2_ex(:, end), xc2);

    err_h_vec2(i) =  1/sqrt(N) * norm(h2(:, end)' - h2ex_interp);
    err_m_vec2(i) =  1/sqrt(N) * norm(m2(:, end)' - m2ex_interp);

end


% Plot the error
figure(4)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec2 , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, 10 * delta_x_vec, "--", "Linewidth", 2)
loglog(delta_x_vec, 100 * delta_x_vec.^2, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^3, "--", "Linewidth", 2)
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=1\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", ...
    "interpreter", "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on


subplot(2,1,2)
loglog(delta_x_vec, err_m_vec2, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, 10 * delta_x_vec, "--", "Linewidth", 2)
loglog(delta_x_vec, 100 * delta_x_vec.^2, "--", "Linewidth", 2)
loglog(delta_x_vec, 10 * delta_x_vec.^3, "--", "Linewidth", 2)
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=1\)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "\(\Delta x^3\)", ...
    "interpreter", "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)