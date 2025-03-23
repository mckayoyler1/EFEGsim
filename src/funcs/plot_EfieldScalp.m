function plot_EfieldScalp(cartPositions, gradients, values, gradientType)
% PLOT_EFIELD Plots gradient vectors at specified points in 3D space.
%
% Inputs:
%   cartPositions - Nx3 matrix of Cartesian coordinates [x, y, z].
%   gradients     - Nx3 matrix of gradients, either in Cartesian or spherical coordinates.
%   values        - Nx1 vector of scalar values at each point (optional, for color coding).
%   gradientType  - String indicating gradient type: 'cartesian' (default) or 'spherical'.
%
% The function visualizes the gradient vectors as arrows, with their magnitudes and
% directions in 3D space. The colors represent the values of the scalar field.

% Validate inputs
if size(cartPositions, 2) ~= 3
    error('Positions must be an Nx3 matrix of Cartesian coordinates [x, y, z].');
end
if size(gradients, 2) ~= 3
    error('Gradients must be an Nx3 matrix.');
end
if nargin < 3 || isempty(values)
    values = zeros(size(cartPositions, 1), 1); % Default to no color coding
end
if nargin < 4
    gradientType = 'cartesian'; % Default to Cartesian gradients
end
mat = [toSphere(cartPositions) gradients values];


gradients = gradients(abs(mat(:,1))<10+1e-3,:);
values = values(abs(mat(:,1))<10+1e-3,:);
cartPositions = cartPositions(abs(mat(:,1))<10+1e-3,:);
mat = mat(abs(mat(:,1))<10+1e-3,:);
% Convert gradients to Cartesian if they are in spherical coordinates
if strcmpi(gradientType, 'spherical')
    % Convert positions to spherical coordinates
    r = mat(:, 1);
    theta = mat(:, 2);
    phi = mat(:, 3);

    % Convert gradients from spherical to Cartesian
    % dr is radial, dtheta is polar, dphi is azimuthal
    gx = gradients(:, 1) .* sin(theta) .* cos(phi) + ...
         (1 ./ r) .* gradients(:, 2) .* cos(theta) .* cos(phi) - ...
         (1 ./ (r .* sin(theta))) .* gradients(:, 3) .* sin(phi);

    gy = gradients(:, 1) .* sin(theta) .* sin(phi) + ...
         (1 ./ r) .* gradients(:, 2) .* cos(theta) .* sin(phi) + ...
         (1 ./ (r .* sin(theta))) .* gradients(:, 3) .* cos(phi);

    gz = gradients(:, 1) .* cos(theta) - ...
         (1 ./ r) .* gradients(:, 2) .* sin(theta);
else
    % Assume gradients are already in Cartesian coordinates
    gx = gradients(:, 1);
    gy = gradients(:, 2);
    gz = gradients(:, 3);
end

% % Normalize gradient vectors for better visualization
% g_magnitude = sqrt(gx.^2 + gy.^2 + gz.^2);
% gx = gx ./ g_magnitude;
% gy = gy ./ g_magnitude;
% gz = gz ./ g_magnitude;


% Extract Cartesian positions
x = cartPositions(:, 1);
y = cartPositions(:, 2);
z = cartPositions(:, 3);

% Plot the scalar field points with color
% ft_plot_topo3d(cartPositions, values, 'neighbourdist', inf);
% hold on;
scatter3(x,y,z, 50, values, 'filled');

% Plot gradient vectors as arrows
quiver3(x,y,z, gx, gy, gz, 2, 'k'); % Black arrows for gradients
colorbar; % Add colorbar for scalar field values
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Gradient Vectors in 3D Space');
axis equal;
grid on;
hold off;
end
