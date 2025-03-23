classdef FieldCalculations
    properties
        Locations    % Nx3 matrix of [x, y, z] positions
        Potentials   % Nx1 vector of potentials at each location
        Interpolant  % Scattered interpolant for the potential field
    end

    methods
        function obj = FieldCalculations(locations, potentials, config) %config gives Interpolation Parameters
            % Constructor to initialize FieldCalculations
            if size(locations, 2) ~= 3
                error('Locations must be an Nx3 matrix of [x, y, z] positions.');
            end

            if length(potentials) ~= size(locations, 1)
                error('Number of potentials must match the number of locations.');
            end
            interpolationMethod = config.interpolationParams.method;
            extrapolationMethod = config.interpolationParams.extrapolation;

            obj.Locations = locations;
            obj.Potentials = potentials;
            obj.Interpolant = scatteredInterpolant(locations(:, 1), locations(:, 2), locations(:, 3), potentials, interpolationMethod, extrapolationMethod);
        end

        function [Fx, Fy, Fz] = computeGradient(obj, newPoints)
            % Compute the gradient of potentials at given locations or new points
            % Input:
            %   newPoints (optional) - Mx3 matrix of [x, y, z] positions to compute gradient
            % Output:
            %   Fx, Fy, Fz - Gradients of the potential field along x, y, and z

            if nargin < 2 || isempty(newPoints)
                newPoints = obj.Locations; % Default to original locations
            end

            if size(newPoints, 2) ~= 3
                error('newPoints must be an Nx3 matrix of [x, y, z] positions.');
            end

            % Compute the gradient at the specified points
            [Fx, Fy, Fz] = obj.Interpolant.gradient(newPoints(:, 1), newPoints(:, 2), newPoints(:, 3));
        end
    end
end