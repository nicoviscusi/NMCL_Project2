clear
close all
clc

%%% Code by Francesco Sala and Nicolo' Viscusi %%%

% Set to true if you want to see the animation of the solutions over time
animation = "True";

%% First set of initial conditions

% Spatial domain
xspan = [0 2];

% Temporal domain
tspan = [0 0.5];

% Initial conditions
h01 = @(x) 1 - 0.1 * sin(pi * x);
m01 = @(x) 0;

% Source function
S1 = @(x, t) [0;
    0];

% We will use periodic boundary condition option ('peri')
bc = 'peri';


% We generate a reference solution
[h1_ex, m1_ex, tvec1_ex, xvec1_ex] = conservative_scheme(xspan, tspan, ...
    3000, 6000, h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);


% We now proceed with a less refined solution
N = 250;

% Number of time steps
CFL = 0.5;
K = N / CFL;

% Solve the problem
[h1, m1, tvec1, xvec1] = conservative_scheme(xspan, tspan, N, K, ...
    h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);


% We visualize the solution
if animation == "True"
    figure(1)
    for i = 1 : 20 : length(tvec1)

        subplot(2, 1, 1)
        plot(xvec1, h1(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec1(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.9 1.1]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec1, m1(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec1(i))], ...
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



%% Error analysis, initial condition 1

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(6:10);

% Note that we cannot solve for small values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec1 = zeros(size(N_vec));
err_m_vec1 = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);

    K = N / CFL;

    [h1, m1, ~, xvec1_err] = conservative_scheme(xspan, tspan, N, K, ...
        h01, m01, @lax_friedrichs_flux, @flux_phys, S1, bc);
    
    % We now want to compare h1(:, end) with h1_ex(:, end), 
    % but this second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xvec1
    h1_interp = interp1(xvec1_err, h1(:, end), xvec1_ex);
    m1_interp = interp1(xvec1_err, m1(:, end), xvec1_ex);
    err_h_vec1(i) = norm(h1_interp' -h1_ex(:, end)); 
    err_m_vec1(i) = norm(m1_interp' - m1_ex(:, end));

end


% Plot error
figure(2)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec1, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=0.5\) (case 1)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on

subplot(2,1,2)
loglog(delta_x_vec, err_m_vec1, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=0.5\) (case 1)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)



%% Second set of initial conditions

% Initial conditions
h02 = @(x) 1 - 0.2 * sin(2 * pi * x);
m02 = @(x) 0.5;

% Source term
S2 = @(x, t) [0;
    0];

% All the other parameters remain the same...

% First, a refined solution as reference "exact"
[h2_ex, m2_ex, tvec2_ex, xvec2_ex] = conservative_scheme(xspan, ...
    tspan, 3000, 6000, h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);

% Solve the problem on a less refined mesh
N = 100;
K = 200;
[h2, m2, tvec2, xvec2] = conservative_scheme(xspan, tspan, N, K, ...
    h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);



% We visualize the solution
if animation == "True"
    figure(3)
    for i = 1 : 20 : length(tvec2)

        subplot(2, 1, 1)
        plot(xvec2, h2(:, i), 'LineWidth', 2)
        title(['$h(x, t)$ at $t = $', num2str(tvec2(i))], ...
            'Interpreter', 'latex')
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$h(x, t)$', 'Interpreter', 'latex')
        grid on
        xlim([0 2]);
        ylim([0.8 1.2]);
        set(gca, 'Fontsize', 20)
        drawnow

        subplot(2, 1, 2)
        plot(xvec2, m2(:, i), 'LineWidth', 2)
        title(['$m(x, t)$ at $t = $', num2str(tvec2(i))], ...
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



%% Error analysis, initial condition 2

% We solve the same problem for different values of \Delta x
delta_x_vec =  2.^-(6:10);

% Note that we cannot solve for small values of delta_x, because we would
% need a too large matrix to store the solutions h and m
N_vec = (xspan(2) - xspan(1)) ./ delta_x_vec ;
err_h_vec2 = zeros(size(N_vec));
err_m_vec2 = zeros(size(N_vec));

for i=1:length(N_vec)
    N = N_vec(i);

    K = N / CFL;

    [h2, m2, ~, xvec2_err] = conservative_scheme(xspan, tspan, N, K, ...
        h02, m02, @lax_friedrichs_flux, @flux_phys, S2, bc);

    % We now want to compare h1(:, end) with h1_ex(:, end), but this 
    % second vector is defined on a different grid xvec1_ex
    % We interpolate h1_ex(:, end) on the grid xvec1
    h2_interp = interp1(xvec2_err, h2(:, end), xvec2_ex);
    m2_interp = interp1(xvec2_err, m2(:, end), xvec2_ex);

    err_h_vec2(i) = norm(h2_interp' - h2_ex(:, end)); 
    err_m_vec2(i) = norm(m2_interp' - m2_ex(:, end));

end


% Plot the error
figure(4)

subplot(2,1,1)
loglog(delta_x_vec, err_h_vec2 , "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(h(x,t)\) at \(t=0.5\) (case 2)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex",  "location", "best")
set(gca, 'Fontsize', 20)
grid on

subplot(2,1,2)
loglog(delta_x_vec, err_m_vec2, "o-", "Linewidth", 2)
hold on
loglog(delta_x_vec, delta_x_vec, "--", delta_x_vec, delta_x_vec.^2, "--")
xlabel('$\Delta x$', 'Interpreter', 'latex')
ylabel("$\|e\|_2$", "Interpreter","latex")
title("Error on \(m(x,t)\) at \(t=0.5\) (case 2)", "Interpreter","latex")
legend("Error", "\(\Delta x\)", "\(\Delta x^2\)", "interpreter", ...
    "latex", "location", "best")
grid on
set(gca, 'Fontsize', 20)