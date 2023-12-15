function f = flux_phys(q)

% FLUX_PHYS - Computes the physical flux function for the shallow
%             water equations.
%
%   f = flux_phys(q)
%
% INPUTS:
%   q   - Vector of state variables [h, m], where
%           h: Water depth
%           m: Discharge
%
% OUTPUT:
%   f   - Vector representing the physical flux corresponding to the input
%         state q.
%
% DESCRIPTION:
%   This function calculates the physical flux for the shallow water
%   equations.
%
% Authors: [Francesco Sala, Nicolo' Viscusi]
% December 2023

g = 1;
h = q(1);
m = q(2);

f = [m;
     m.^2./h + 1/2*g*h.^2];

return