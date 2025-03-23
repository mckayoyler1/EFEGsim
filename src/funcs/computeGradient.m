function gradients = computeGradient(points, scalarField, config, newPoints)
% computeGradient computes the gradient of a scalar field at given points.
% 
% Inputs:
%   points      - Nx3 array of points (x, y, z).
%   scalarField - Nx1 array of scalar values corresponding to the points.
%   newPoints   - Mx3 array of new points where the gradient should be evaluated (optional).
% 
% Output:
%   gradients   - Nx3 or Mx3 array of gradient vectors at the specified points.

% Check inputs
if nargin < 2
    error('At least two inputs (points and scalarField) are required.');
end

if size(points, 2) ~= 3
    error('Points must be an Nx3 array of (x, y, z) coordinates.');
end

if size(scalarField, 1) ~= size(points, 1)
    error('Scalar field must have the same number of rows as points.');
end

if nargin < 4
    newPoints = points; % Default to calculating gradients at input points
end

interpolationMethod = config.interpolationParams.method;
extrapolationMethod = config.interpolationParams.extrapolation;
fprintf("Computing gradient -- interpolation method: '%s' | extrapolation method: '%s'.\n", interpolationMethod, extrapolationMethod);
% Compute the gradient using scattered interpolation
V = scatteredInterpolant(points, scalarField, interpolationMethod, extrapolationMethod);

% Compute partial derivatives in x, y, z
[Fx, Fy, Fz] = gradientComponents(V, newPoints);
gradients = [Fx, Fy, Fz];

end

function [Fx, Fy, Fz] = gradientComponents(F, points)
% Compute partial derivatives in x, y, and z directions
    delta = 1e-5; % Small increment for numerical differentiation
    
    Fx = (F(points + [delta, 0, 0]) - F(points - [delta, 0, 0])) / (2 * delta);
    Fy = (F(points + [0, delta, 0]) - F(points - [0, delta, 0])) / (2 * delta);
    Fz = (F(points + [0, 0, delta]) - F(points - [0, 0, delta])) / (2 * delta);
end
