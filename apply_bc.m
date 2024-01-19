function U = apply_bc(Ui, bc, m)
% APPLY_BC - Apply boundary conditions to a given solution matrix.
%
%   U = APPLY_BC(Ui, bc, m) applies boundary conditions to the input matrix
%   UI based on the specified boundary condition type 'BC' and the number of
%   ghost points 'M'. The function returns the modified matrix 'U' with
%   added ghost points.
%
%   Input:
%       Ui   - Input matrix without boundary conditions.
%       bc   - String specifying the type of boundary condition. 
%              Supported values: 'peri' (periodic), 'open' (open boundary).
%       m    - Number of ghost points to be added on each side.
%
%   Output:
%       U    - Matrix with applied boundary conditions.
%
% Function available on Moodle for MATH-459 course
% January 2024


switch bc

    case 'peri'

        U = [Ui(:, end-m+1 : end), Ui, Ui(:,1:m)];

    case 'open'

        U = [repmat(Ui(:,1), 1, m), Ui, repmat(Ui(:,end),1,m)];

end

end