function plotSymbolicEField(EField, points)
    % plotSymbolicEField Plots the symbolic electric field as vectors at specified points
    %
    % Inputs:
    %   EField - 1x3 cell array of symbolic expressions for [Ex, Ey, Ez]
    %   points - Nx3 matrix of [x, y, z] points where the E-field is evaluated

    % Validate inputs
    if size(points, 2) ~= 3
        error('Points must be an Nx3 matrix of [x, y, z] coordinates.');
    end

[Ex_values, Ey_values, Ez_values] = deal(evaluateSymE(EField, points));

    % Plot the electric field vectors
    figure;
    quiver3(points(:, 1), points(:, 2), points(:, 3), ...
            Ex_values, Ey_values, Ez_values, 10);
    hold on;
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Electric Field Vector Field');
    grid on;
    axis equal;
end
