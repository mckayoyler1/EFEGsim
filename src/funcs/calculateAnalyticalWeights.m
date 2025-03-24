function data = calculateAnalyticalWeights(data)
% CALCULATEANALYTICALWEIGHTS Calculates analytical weights for points on a uniform spherical grid
%
% Inputs:
%   data - Data structure containing:
%          .sphpos - Spherical coordinates
%          .pos - Cartesian coordinates
%          .grid_type - Must be 'uniform' for this function
%          .cfg.sensorArray.numTheta - Number of theta points
%          .cfg.sensorArray.numPhi - Number of phi points
%          .cfg.headshape.radius - Radius of the spherical shell
%
% Outputs:
%   data - Updated data structure with added field:
%          .voronoi_weights - Vector of weights for numerical integration
%
% Note: This function can ONLY be used with a uniform grid in spherical coordinates.
%       It will throw an error if used with 'random' or 'modulated' grids.

% % Check if grid type is 'uniform'
% if ~isfield(data, 'grid_type')
%     data.grid_type = 'uniform'; % Default to uniform if not specified
% elseif ~strcmpi(data.grid_type, 'uniform') || ~strcmpi(data.grid_type, 'modulated')
%     error('calculateAnalyticalWeights can only be used with uniform grids (grid_type = ''uniform''). Current grid_type is ''%s''.', data.grid_type);
% end

% Get spherical coordinates
sphpos = data.sphpos;
theta_values = sphpos(:,2);
% phi_values = spherical_coords(:,3);
r = data.cfg.headshape.radius;
% Get grid dimensions
n_theta = data.cfg.sensorArray.numTheta;
n_phi = data.cfg.sensorArray.numPhi;

% Calculate spacing
dtheta = max(theta_values) / (n_theta - 1);
dphi = 2*pi / (n_phi - 1);

% Initialize weights array
weights = zeros(size(theta_values));

% Get unique theta values
unique_thetas = unique(theta_values);

% For each unique theta, assign the appropriate weight to all points with that theta
for i = 1:length(unique_thetas)
    current_theta = unique_thetas(i);
    
    % Find all points with this theta
    theta_mask = abs(theta_values - current_theta) < 1e-10;
    
    % Calculate the weight for this theta
    theta_weight = r^2 * sin(current_theta) * dtheta * dphi;
    
    % Assign this weight to all points with this theta
    weights(theta_mask) = theta_weight;
end

% Normalize to hemisphere surface area

% weights = 4 * pi * r * weights / sum(weights);

% Store the weights
data.weights = weights;

% ----- Special handling for pole points -----
% Identify pole points (where θ is very close to 0)
is_pole = abs(sphpos(:,2)) < 1e-10;
num_pole_points = sum(is_pole);

if num_pole_points > 1
    fprintf('Found %d pole points, adjusting weights\n', num_pole_points);
    
    % Find the indices of pole points
    pole_indices = find(is_pole);
    
    % Keep one pole point and adjust its weight
    first_pole = pole_indices(1);
    other_poles = pole_indices(2:end);
    
    % Sum the weights of all pole points
    total_pole_weight = sum(data.weights(pole_indices));
    
    % Assign the total weight to the first pole point
    data.weights(first_pole) = total_pole_weight;
    
    % Remove the redundant pole points
    data.sphpos(other_poles,:) = [];
    data.pos(other_poles,:) = [];
    data.weights(other_poles) = [];
    
    % Update the number of points
    data.num_points = size(data.pos, 1);
    
    fprintf('After pole handling: %d points remain\n', data.num_points);
end

end