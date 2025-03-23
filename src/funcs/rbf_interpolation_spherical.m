function interpolated_values = rbf_interpolation_spherical(points, values, query_points, kernel, epsilon)
% RBF_INTERPOLATION Interpolates values using Radial Basis Function (RBF) interpolation.
%
% Inputs:
%   points         - NxD matrix of input points (N points in D-dimensional space).
%   values         - Nx1 vector of scalar values at each input point.
%   query_points   - MxD matrix of query points (M points in D-dimensional space).
%   kernel         - RBF kernel type: 'gaussian', 'multiquadric', 
%                    'inverse_multiquadric', or 'thin_plate'.
%   epsilon        - Shape parameter for the RBF (positive scalar, required for 
%                    most kernels except 'thin_plate').
%
% Outputs:
%   interpolated_values - Mx1 vector of interpolated values at the query points.

% Validate inputs
if nargin < 4
    error('At least four arguments are required: points, values, query_points, and kernel.');
end
if nargin < 5 && ~strcmp(kernel, 'thin_plate')
    error('The shape parameter epsilon is required for the selected kernel.');
end

% Number of input points
N = size(points, 1);

% Distance function
distance = @(x1, x2) sqrt(sum((x1 - x2).^2, 2)); % Euclidean distance

% Define RBF kernels
switch kernel
    case 'gaussian'
        rbf = @(r) exp(-(epsilon * r).^2);
    case 'multiquadric'
        rbf = @(r) sqrt(1 + (epsilon * r).^2);
    case 'inverse_multiquadric'
        rbf = @(r) 1 ./ sqrt(1 + (epsilon * r).^2);
    case 'thin_plate'
        rbf = @(r) r.^2 .* log(r + (r == 0)); % Handle r == 0 to avoid log(0)
    otherwise
        error('Unsupported RBF kernel. Choose: gaussian, multiquadric, inverse_multiquadric, or thin_plate.');
end

% Compute the pairwise distance matrix between input points
distance_matrix = zeros(N, N);
for i = 1:N
    distance_matrix(i, :) = distance(points(i, :), points);
end

% Compute the RBF matrix
rbf_matrix = rbf(distance_matrix);

% Solve for the weights
weights = rbf_matrix \ values;

% Number of query points
M = size(query_points, 1);

% Interpolate values at query points
interpolated_values = zeros(M, 1);
for i = 1:M
    % Compute distances between the query point and all input points
    query_distances = distance(query_points(i, :), points);
    
    % Compute the RBF values
    query_rbf = rbf(query_distances);
    
    % Compute the interpolated value
    interpolated_values(i) = sum(weights .* query_rbf);
end
end
