classdef DipoleDistribution < handle
    properties
        Dipoles            % Array of Dipole objects
        Locations          % Nx3 matrix of dipole locations [x, y, z]
        Moments            % Nx3 matrix of dipole moments [px, py, pz]
        SymbolicPotential  % Symbolic expression for total potential
    end

    methods
        function obj = DipoleDistribution(dipoles)
            % Constructor: Initialize with a list of Dipole objects
            if nargin > 0 && ~isempty(dipoles)
                obj.Dipoles = dipoles;
                obj.updateLocationsAndMoments();
            else
                obj.Dipoles = DipoleDistribution.empty; % Empty array of Dipole objects
                obj.Locations = [];
                obj.Moments = [];
            end
        end

        function addDipole(obj, dipole)
            % Add a single Dipole object to the distribution
            if ~isa(dipole, 'Dipole')
                error('Input must be an instance of the Dipole class.');
            end
            obj.Dipoles = [obj.Dipoles, dipole];
            obj.updateLocationsAndMoments();
        end

        function removeDipole(obj, index)
            % Remove a dipole by its index
            if index < 1 || index > numel(obj.Dipoles)
                error('Index out of bounds.');
            end
            obj.Dipoles(index) = [];
            obj.updateLocationsAndMoments();
        end

        function updateLocationsAndMoments(obj)
            % Update Locations and Moments from the list of Dipoles
            numDipoles = numel(obj.Dipoles);
            obj.Locations = zeros(numDipoles, 3);
            obj.Moments = zeros(numDipoles, 3);

            for i = 1:numDipoles
                obj.Locations(i, :) = obj.Dipoles(i).Position;
                obj.Moments(i, :) = obj.Dipoles(i).Moment;
            end
        end


        function visualize(obj)
            % Visualize the dipole distribution
            figure;
            hold on;

            for i = 1:numel(obj.Dipoles)
                dipole = obj.Dipoles(i);
                % Plot dipole location as a point
                scatter3(dipole.Position(1), dipole.Position(2), dipole.Position(3), 100, 'filled');
                % Plot dipole moment as a vector
                quiver3(dipole.Position(1), dipole.Position(2), dipole.Position(3), ...
                        dipole.Moment(1), dipole.Moment(2), dipole.Moment(3), ...
                        'LineWidth', 2, 'MaxHeadSize', 1);
            end

            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Dipole Distribution');
            axis equal;
            grid on;
            hold off;
        end
    end
end
