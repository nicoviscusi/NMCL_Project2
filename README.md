# Numerical Methods for Conservation Laws: Project 2

Welcome to the repository for the 2nd project of EPFL MATH-459: Numerical Methods for Conservation Laws. This project is dedicated to the implementation of WENO schemes for a specific Partial Differential Equation (PDE) problem, that is the shallow water equations, as outlined in the assignment. Note that this project exploits a few functions implemented by Jan S. Hesthaven.

## Code Structure

The project code is organized into two main files:

1. `part2_3.m`
2. `part2_4.m`

Each of these files executes a numerical simulation for the PDE problem specified in the project assignment. The user can control whether to visualize an animation of the solution by modifying the first line of each file, setting `animation` to either "true" or "false."

## Core Functions

### 1. `solver.m`

The generic numerical conservative scheme is implemented in this function. It takes several inputs, including a physical flux and a numerical flux, and plays a central role in the numerical simulations conducted in the project.

### 2. `evalRHS.m`

The RHS of the hyperbolic equation is evaluated at each timestep, and the reconstrucion of the cells values is performed thanks to the WENO scheme.

### 3.  `WENO.m`

WENO reconstruction is performed by this Matlab script.


## Usage

To execute a specific part of the project, run the corresponding `part2_X.m` file. The code can run for any value of k, i.e. the order of accuracy of the WENO scheme. Adjust the `animation` variable in the first lines to control the visualization of the solution. The code will perform an error analysis of the WENO scheme, using as reference solution the exact solution if available, or a solution computed on a very fine mesh.

## Project Contributors

This project was realized by:

- Francesco Sala
- Nicolo' Viscusi

Completed in January 2024.

If you have any questions or suggestions, please don't hesitate to reach out to the contributors.

