function data = getShellPoints(data)
% GETSHELLPOINTS Generates points on a spherical shell
%
% Inputs:
%   data - Data structure with fields:
%       .cfg.headshape.radius - Radius of the spherical shell
%       .cfg.sensorArray.thetaRange - [min_theta, max_theta] in radians
%       .cfg.sensorArray.phiRange - [min_phi, max_phi] in radians
%       .cfg.sensorArray.numTheta - Number of theta points
%       .cfg.sensorArray.numPhi - Number of phi points
%       .grid_type - String specifying grid type:
%           'uniform' - Regular grid (default)
%           'random' - Random points on sphere
%           'modulated' - Modulated radius (±5%)
%           'all' - Generate all three grid types
%       .mod_percent - Percentage to modulate radius (default: 5)
%
% Outputs:
%   data - Updated data structure with:
%       .pos - Nx3 array of cartesian [x, y, z] coordinates for the selected grid type
%       .sphpos - Nx3 array of spherical [r, theta, phi] coordinates for the selected grid type
%       .grid - Struct containing the selected grid information
%       Additionally, if grid_type='all', the following fields are also added:
%           .uniform_grid - Struct with uniform grid
%           .random_grid - Struct with random grid
%           .modulated_grid - Struct with modulated radius grid

% Set default grid type if not provided
if ~isfield(data, 'grid_type')
    data.grid_type = 'uniform';
end

if ~isfield(data, 'mod_percent')
    data.mod_percent = 5;
end

% Extract parameters from data struct
r = data.cfg.headshape.radius;
theta_min = data.cfg.sensorArray.thetaRange(1);
theta_max = data.cfg.sensorArray.thetaRange(2);
phi_min = data.cfg.sensorArray.phiRange(1);
phi_max = data.cfg.sensorArray.phiRange(2);

% Calculate number of points in grid
n_theta = data.cfg.sensorArray.numTheta;
n_phi = data.cfg.sensorArray.numPhi;
total_points = n_theta * n_phi;

[theta, phi] = ndgrid(linspace(theta_min, theta_max, n_theta), ...
                      linspace(phi_min, phi_max, n_phi));

% Reshape to column vectors
theta = theta(:);
phi = phi(:);

% Generate the grid(s) based on grid_type
grid_type = lower(data.grid_type);

% Check if we need a uniform grid
if strcmp(grid_type, 'uniform') || strcmp(grid_type, 'all')
    uniform_grid = struct();
    uniform_grid.num_points = total_points;
    
    % Create spherical coordinates for uniform grid
    uniform_grid.sphpos = [repmat(r, total_points, 1), theta, phi];
    
    % Convert to Cartesian coordinates
    uniform_grid.pos = spherical_to_cartesian(uniform_grid.sphpos);
    
    % If generating only uniform grid or 'all' is specified
    if strcmp(grid_type, 'uniform')
        % Add to data struct as main position references
        data.pos = uniform_grid.pos;
        data.sphpos = uniform_grid.sphpos;
        data.grid = uniform_grid;
    end
    
    % If we're generating all types, save it as uniform_grid
    if strcmp(grid_type, 'all')
        data.uniform_grid = uniform_grid;
    end
end

% Check if we need a random grid
if strcmp(grid_type, 'random') || strcmp(grid_type, 'all')
    random_grid = struct();
    random_grid.num_points = total_points;
    
    % Generate random points with uniform distribution on sphere
    spherical_coords = zeros(total_points, 3);
    
    % All points have the same radius
    spherical_coords(:,1) = r;
    
    % Generate random angles with correct distribution for uniform coverage
    % For uniform distribution on sphere, theta should be acos(2*rand-1)
    spherical_coords(:,2) = acos(1 - 2*rand(total_points, 1)); % theta in [0, pi]
    spherical_coords(:,3) = 2*pi*rand(total_points, 1);        % phi in [0, 2*pi]
    
    % Adjust angles to specified ranges if needed
    if theta_min > 0 || theta_max < pi
        theta_range = theta_max - theta_min;
        spherical_coords(:,2) = theta_min + theta_range * rand(total_points, 1);
    end
    
    if phi_min > 0 || phi_max < 2*pi
        phi_range = phi_max - phi_min;
        spherical_coords(:,3) = phi_min + phi_range * rand(total_points, 1);
    end
    
    random_grid.sphpos = spherical_coords;
    random_grid.pos = spherical_to_cartesian(random_grid.sphpos);
    
    % If generating only random grid
    if strcmp(grid_type, 'random')
        % Add to data struct as main position references
        data.pos = random_grid.pos;
        data.sphpos = random_grid.sphpos;
        data.grid = random_grid;
    end
    
    % If we're generating all types, save it as random_grid
    if strcmp(grid_type, 'all')
        data.random_grid = random_grid;
    end
end

% Check if we need a modulated radius grid
if strcmp(grid_type, 'modulated') || strcmp(grid_type, 'all')
    modulated_grid = struct();
    modulated_grid.num_points = total_points;
    
    % Create base spherical coordinates (same angular coordinates as uniform grid)
    modulated_spherical = [repmat(r, total_points, 1), theta, phi];
    
    % Calculate modulation amount
    mod_range = data.mod_percent / 100 * r;
    
    % Create random modulation for each point
    modulations = 2 * mod_range * rand(total_points, 1) - mod_range;
    
    % Apply modulation to radius
    modulated_spherical(:,1) = r + modulations;
    
    modulated_grid.sphpos = modulated_spherical;
    modulated_grid.pos = spherical_to_cartesian(modulated_grid.sphpos);
    
    % If generating only modulated grid
    if strcmp(grid_type, 'modulated')
        % Add to data struct as main position references
        data.pos = modulated_grid.pos;
        data.sphpos = modulated_grid.sphpos;
        data.grid = modulated_grid;
    end
    
    % If we're generating all types, save it as modulated_grid
    if strcmp(grid_type, 'all')
        data.modulated_grid = modulated_grid;
    end
end

% If grid_type is 'all', we need to pick one for the main position references
% Let's default to uniform for the main reference
if strcmp(grid_type, 'all')
    data.pos = data.uniform_grid.pos;
    data.sphpos = data.uniform_grid.sphpos;
    data.grid = data.uniform_grid;
end

% Also store the number of points for convenience
data.num_points = total_points;

end

% Helper function to convert spherical to Cartesian coordinates
function cartesian = spherical_to_cartesian(spherical)
    r_vals = spherical(:,1);
    theta_vals = spherical(:,2);
    phi_vals = spherical(:,3);
    
    cartesian = zeros(size(spherical));
    cartesian(:,1) = r_vals .* sin(theta_vals) .* cos(phi_vals); % x
    cartesian(:,2) = r_vals .* sin(theta_vals) .* sin(phi_vals); % y
    cartesian(:,3) = r_vals .* cos(theta_vals);                  % z
end