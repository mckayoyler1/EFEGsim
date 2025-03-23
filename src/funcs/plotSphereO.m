function plotSphereO(radius, alpha)
% Parameters
center = [0, 0, 0];  % Center of the sphere (x, y, z)

% Create sphere data
[theta, phi] = meshgrid(linspace(0, pi, 50), linspace(0, 2*pi, 100));
X = radius * sin(theta) .* cos(phi) + center(1);  % X-coordinates
Y = radius * sin(theta) .* sin(phi) + center(2);  % Y-coordinates
Z = radius * cos(theta) + center(3);  % Z-coordinates

% Plot the sphere
figure;
sphere_surface = surf(X, Y, Z);

% Set appearance properties
sphere_surface.FaceAlpha = alpha;  % Set transparency
sphere_surface.EdgeColor = 'none';        % Remove mesh edges
sphere_surface.FaceColor = [0.2, 0.6, 0.8];  % Set sphere color (e.g., light blue)

% Adjust axes for proper visualization
axis equal;  % Equal aspect ratio
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Semi-Transparent Sphere');
grid on;
view(3);  % Set 3D view

end