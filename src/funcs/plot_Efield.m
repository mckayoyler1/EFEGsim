function plot_Efield(cartPositions, vectors, gradientType)
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
if size(vectors, 2) ~= 3
    error('Gradients must be an Nx3 matrix.');
end
if nargin < 5 || isempty(values)
    values = zeros(size(cartPositions, 1), 1); % Default to no color coding
end
if nargin < 3
    gradientType = 'cartesian'; % Default to Cartesian gradients
end

% Convert gradients to Cartesian if they are in spherical coordinates
if strcmpi(gradientType, 'spherical')
    % Convert positions to spherical coordinates
    sphPositions = toSphere(cartPositions);
    r = sphPositions(:, 1);
    theta = sphPositions(:, 2);
    phi = sphPositions(:, 3);

    % Convert gradients from spherical to Cartesian
    % dr is radial, dtheta is polar, dphi is azimuthal
    gx = vectors(:, 1) .* sin(theta) .* cos(phi) + ...
         (1 ./ r) .* vectors(:, 2) .* cos(theta) .* cos(phi) - ...
         (1 ./ (r .* sin(theta))) .* vectors(:, 3) .* sin(phi);

    gy = vectors(:, 1) .* sin(theta) .* sin(phi) + ...
         (1 ./ r) .* vectors(:, 2) .* cos(theta) .* sin(phi) + ...
         (1 ./ (r .* sin(theta))) .* vectors(:, 3) .* cos(phi);

    gz = vectors(:, 1) .* cos(theta) - ...
         (1 ./ r) .* vectors(:, 2) .* sin(theta);
else
    % Assume gradients are already in Cartesian coordinates
    gx = vectors(:, 1);
    gy = vectors(:, 2);
    gz = vectors(:, 3);
end


% Extract Cartesian positions
x = cartPositions(:, 1);
y = cartPositions(:, 2);
z = cartPositions(:, 3);

% % Normalize gradient vectors for better visualization
% g_magnitude = sqrt(gx.^2 + gy.^2 + gz.^2);
% gx = gx ./ g_magnitude;
% gy = gy ./ g_magnitude;
% gz = gz ./ g_magnitude;


% % Plot the scalar field points with color
% ft_plot_topo3d(cartPositions, values, 'neighbourdist', inf);
% hold on;

% Plot gradient vectors as arrows
quiver3(x, y, z, gx, gy, gz, 2, 'k'); % Black arrows for gradients
colorbar; % Add colorbar for scalar field values
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Gradient Vectors in 3D Space');
axis equal;
grid on;
hold off;
end
