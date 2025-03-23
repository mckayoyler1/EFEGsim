function plotDipoles(dipoles)
    % visualizeDipoles Plots dipole locations and moments in 3D space.
    % Inputs:
    %   dipoles - Array of Dipole objects, each with Position and Moment properties.

    % Check if the input is valid
    if isempty(dipoles)
        error('The input list of dipoles is empty.');
    end

    % Initialize figure
    figure;
    hold on;
    grid on;
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(3);
    title('Dipole Locations and Moments');

    % Plot each dipole
    for i = 1:numel(dipoles)
        % Extract position and moment
        pos = dipoles(i).Position;  % [x, y, z]
        moment = .1 * dipoles(i).Moment/norm(dipoles(i).Moment); % [mx, my, mz]

        % Plot position as a red dot
        scatter3(pos(1), pos(2), pos(3), 50, 'r', 'filled');

        % Plot moment as a blue vector
        quiver3(pos(1), pos(2), pos(3), moment(1), moment(2), moment(3), 0, ...
            'Color', 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    end

    legend({'Dipole Position', 'Dipole Moment'}, 'Location', 'best');
    hold off;
end