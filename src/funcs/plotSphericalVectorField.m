function plotSphericalVectorField(sphPoints, sphVectors, scaleFactor)
% plotSphericalVectorField Plots a vector field given in spherical coordinates.
%
%   plotSphericalVectorField(sphPoints, sphVectors) plots the vector field
%   where:
%       - sphPoints is an M-by-3 matrix with rows [r, theta, phi] (angles in radians)
%       - sphVectors is an M-by-3 matrix with rows [V_r, V_theta, V_phi]
%
%   plotSphericalVectorField(..., scaleFactor) scales the vector lengths by
%   the provided scaleFactor (default is 1).
%
%   The spherical vector field is converted to Cartesian coordinates using:
%
%       [x, y, z] = [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)]
%
%   and the vector components are transformed by:
%
%       Vx = V_r*sin(theta)*cos(phi) + V_theta*cos(theta)*cos(phi) - V_phi*sin(phi)
%       Vy = V_r*sin(theta)*sin(phi) + V_theta*cos(theta)*sin(phi) + V_phi*cos(phi)
%       Vz = V_r*cos(theta) - V_theta*sin(theta)
%
%   Example:
%       % Generate sample points and a simple radial vector field.
%       % Points on a sphere of radius 1:
%       [phi, theta] = meshgrid(linspace(0,2*pi,20), linspace(0,pi,10));
%       r = ones(size(theta));
%       sphPoints = [r(:), theta(:), phi(:)];
%
%       % Let the vector field be purely radial:
%       V_r = ones(size(r(:)));
%       V_theta = zeros(size(r(:)));
%       V_phi = zeros(size(r(:)));
%       sphVectors = [V_r, V_theta, V_phi];
%
%       % Plot the vector field.
%       plotSphericalVectorField(sphPoints, sphVectors);
%
%   See also: quiver3

    if nargin < 3 || isempty(scaleFactor)
        scaleFactor = 1;
    end

    % Extract spherical coordinates from sphPoints.
    r = sphPoints(:,1);
    theta = sphPoints(:,2);
    phi = sphPoints(:,3);
    
    % Convert spherical points to Cartesian coordinates.
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
    
    % Extract spherical vector components.
    Vr = sphVectors(:,1);
    Vtheta = sphVectors(:,2);
    Vphi = sphVectors(:,3);
    
    % Convert spherical vector components to Cartesian components.
    Vx = Vr .* sin(theta) .* cos(phi) + Vtheta .* cos(theta) .* cos(phi) - Vphi .* sin(phi);
    Vy = Vr .* sin(theta) .* sin(phi) + Vtheta .* cos(theta) .* sin(phi) + Vphi .* cos(phi);
    Vz = Vr .* cos(theta) - Vtheta .* sin(theta);
    
    % Apply scaling factor if needed.
    Vx = scaleFactor * Vx;
    Vy = scaleFactor * Vy;
    Vz = scaleFactor * Vz;
    
% Create a figure with white background and a suitable size

% Plot using quiver3 with no automatic scaling
q = quiver3(x, y, z, Vx, Vy, Vz, 0, 'LineWidth', 1.5);  % If you want to plot more stuff afterward

% Label axes; optionally use LaTeX interpreter
xlabel('X','FontSize',14,'Interpreter','tex');
ylabel('Y','FontSize',14,'Interpreter','tex');
zlabel('Z','FontSize',14,'Interpreter','tex');

% Make axes equal for proper 3D proportions
axis equal;

% Optionally set limits (comment out if you want auto-limits)
% xlim([-2 2]);
% ylim([-2 2]);
% zlim([-2 2]);

% Show grid
grid on;

% Put a box around the axes
box on;

% Increase the font size of the axes and tick marks
set(gca, 'FontSize',12, 'LineWidth',1.2);

% Adjust the view if desired (default is 3D, but you can set angles)
view(45,30);

% Optionally, if your data is large or small, you can tweak quiver scaling:
% q.AutoScale = 'on';   % or 'off'
% q.AutoScaleFactor = 0.5;  % for example

% Done

end
