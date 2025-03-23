function interpolatedValues = interpolateValues(points, values, newPoints, interpolation)
% INTERPOLATEVALUES Interpolates the values of a function at given points
%                   to a new set of points. sensorArray -> grid
% Input:
%   points        - Nx3 matrix of original points in Cartesian coordinates (x, y, z).
%   values        - Nx1 vector of function values corresponding to the original points.
%   newPoints     - Mx3 matrix of new points in Cartesian coordinates (x, y, z) where interpolation is desired.
% Output:
%   interpolatedValues - Mx1 vector of interpolated function values at the new points.

% Validate inputs
if size(points, 2) ~= 3
    error('Input "points" must be an Nx3 matrix of Cartesian coordinates.');
end
if size(newPoints, 2) ~= 3
    error('Input "newPoints" must be an Mx3 matrix of Cartesian coordinates.');
end
if size(points, 1) ~= length(values)
    error('The number of rows in "points" must match the length of "values".');
end

switch 
% Use scatteredInterpolant for interpolation
F = scatteredInterpolant(points(:, 1), points(:, 2), points(:, 3), values, 'linear', 'none');

% Evaluate the function at the new points
interpolatedValues = F(newPoints(:, 1), newPoints(:, 2), newPoints(:, 3));
end