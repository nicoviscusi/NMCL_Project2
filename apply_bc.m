% Function returns an extended vector, based on
% the type of boundary condition requested
% Vector is extended by m cells on each side

function U = apply_bc(Ui, bc, m)

comp1 = Ui(1, :);
comp2 = Ui(2, :);

switch bc

    case 'peri'

        U(1, :) = [comp1(end-m+1 : end), comp1, comp1(1:m)];
        U(2, :) = [comp2(end-m+1 : end), comp2, comp2(1:m)];
        %U = [Ui(end-m+1 : end), Ui, Ui(1:m)];

    case 'open'

        U = [repmat(Ui(1), 1, m), Ui, repmat(Ui(end),1,m)];

end

end