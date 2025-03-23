function plotCartesianVectorField(cartPoints, cartVectors, scaleFactor)
% plotCartesianVectorField Plots a vector field given in Cartesian coordinates.
%
%   plotCartesianVectorField(cartPoints, cartVectors) plots the vector field
%   using the points and vector components provided.
%
%   plotCartesianVectorField(..., scaleFactor) scales the arrow lengths by the
%   specified scaleFactor (default is 1).
%
%   INPUTS:
%       cartPoints  - An Mx3 matrix of Cartesian coordinates [x, y, z]
%       cartVectors - An Mx3 matrix of Cartesian vector components [V_x, V_y, V_z]
%       scaleFactor - (Optional) A scalar to scale the arrow lengths (default = 1)
%
%   Example:
%       % Suppose cartPoints and cartVectors are obtained from convertSphericalToCartesian:
%       plotCartesianVectorField(cartPoints, cartVectors, 1);
%
%   This function uses quiver3 to display the vector field.

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 1;
end

% Extract coordinates.
x = cartPoints(:,1);
y = cartPoints(:,2);
z = cartPoints(:,3);

% Extract vector components.
Vx = cartVectors(:,1);
Vy = cartVectors(:,2);
Vz = cartVectors(:,3);

% Plot using quiver3 with no automatic scaling
q = quiver3(x, y, z, Vx, Vy, Vz, 0, 'LineWidth', 1.5);
hold on;  % If you want to plot more stuff afterward

% Label axes; optionally use LaTeX interpreter
xlabel('X','FontSize',14,'Interpreter','tex');
ylabel('Y','FontSize',14,'Interpreter','tex');
zlabel('Z','FontSize',14,'Interpreter','tex');

% Give the plot a title
title('Vector Field','FontSize',16,'Interpreter','tex');

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

%Optionally, if your data is large or small, you can tweak quiver scaling:
%q.AutoScale = 'on';   % or 'off'
% q.AutoScaleFactor = 0.5;  % for example

% Done

end