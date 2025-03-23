function EField = evaluateSymbE(pos, symbE)
%EVALUATESYME Evaluates a symbolic electric field at specific 3D points.
%
%   EField = evaluateSymE(symE, pos)
%
%   Inputs:
%     symE - a cell array of length 3, containing symbolic expressions for the
%            x, y, and z components of an electric field, e.g. {Ex_sym, Ey_sym, Ez_sym}.
%     pos  - an N-by-3 matrix of (x, y, z) coordinates at which to evaluate
%            the field. Each row is a point in 3D space.
%
%   Output:
%     EField - an N-by-3 numeric array containing the evaluated electric field
%              at each row of 'pos'. Specifically, EField(i,:) is [Ex, Ey, Ez]
%              evaluated at pos(i,:).

% Define symbolic variables (x, y, z) for converting symbolic expressions
syms x y z

% Convert each symbolic component (Ex, Ey, Ez) into a numeric MATLAB function
% that accepts x, y, z as inputs.
Ex_func = matlabFunction(symbE{1}, 'Vars', [x, y, z]);
Ey_func = matlabFunction(symbE{2}, 'Vars', [x, y, z]);
Ez_func = matlabFunction(symbE{3}, 'Vars', [x, y, z]);

% Evaluate each component of the electric field at the specified array of points.
% 'pos(:,1), pos(:,2), pos(:,3)' are the x, y, z coordinates for all N points.
Ex = Ex_func(pos(:, 1), pos(:, 2), pos(:, 3));
Ey = Ey_func(pos(:, 1), pos(:, 2), pos(:, 3));
Ez = Ez_func(pos(:, 1), pos(:, 2), pos(:, 3));

% Combine the evaluated components into a single N-by-3 numeric array
EField = [Ex, Ey, Ez];

end
