classdef SensorArray < handle
    properties
        Sensors     % List of Sensor objects
        Locations   % Nx3 matrix of [x, y, z] positions
        MeasuredV  % Nx1 vector of potentials at each location
        SymbolicV
        MeasuredE    % Nx3 matrix of electric field components [Fx, Fy, Fz]
        NumSensors
        Arrangement
    end

    methods
        function obj = SensorArray(config, locations)
            % Constructor to initialize sensor locations
            if nargin < 2 || isempty(locations)
                % No locations provided; determine based on arrangement
                radius = config.headshape.radius;
                thetaRange = config.sensorArray.thetaRange;
                phiRange = config.sensorArray.phiRange;
                numTheta = config.sensorArray.numTheta;
                numPhi = config.sensorArray.numPhi;
                seed = config.seed;

                arrangement = lower(config.sensorArray.arrangement);
                fprintf('Generating sensor locations with arrangement: %s\n', arrangement);
                % Generate locations based on arrangement
                switch arrangement
                    case "grid"
                        if config.sensorArray.multiR == true
                            numSensors = numPhi * numTheta * config.headshape.numR;
                            locations = SensorArray.generateGrid(radius, thetaRange, phiRange, numTheta, numPhi, config.headshape.rRange, config.headshape.numR);
                        else
                            numSensors = numPhi * numTheta;
                            locations = SensorArray.generateGrid(radius, thetaRange, phiRange, numTheta, numPhi);
                        end
                    case "random"
                        numSensors = config.sensorArray.numSensors;
                        locations = SensorArray.generateRandom(radius, thetaRange, phiRange, numSensors, seed);
                    otherwise
                        error("ERROR: Invalid array arrangement '%s'. Must be 'grid' or 'random'.", arrangement);
                end
            else
                % Locations provided; validate input
                if size(locations, 2) ~= 3
                    error('Locations must be an Nx3 matrix of [x, y, z] positions.');
                end
            end

            % Assign locations to the object
            obj.Locations = locations;
            obj.NumSensors = numSensors;
            obj.Arrangement = arrangement;
            obj.Sensors = Sensor.generateSensors(obj);
            obj.MeasuredV = zeros(size(locations, 1), 1); % Initialize potentials
            obj.MeasuredE = zeros(size(locations)); % Initialize EField
            obj.SymbolicV = 0;
        end
        function measurePotential(obj, dipoles)
            % Compute potentials at each sensor location
            % Inputs:
            %   dipoles - Array of Dipole objects
            %   k       - Constant for potential calculation

            potential = zeros(obj.NumSensors, 1);

            for i = 1:numel(dipoles)
                dipole = dipoles(i);
                potential = potential + dipole.V(obj.Locations);
            end
            % Update potentials property and Sensors
            obj.MeasuredV = potential;
            for i = 1:obj.NumSensors
                obj.Sensors(i).Potential = potential(i);
            end
        end
        function computeE(obj)
               V = obj.MeasuredV;



            obj.MeasuredE = Efield;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotSensors(obj)
            % Visualize sensor locations
            scatter3(obj.Locations(:, 1), obj.Locations(:, 2), obj.Locations(:, 3), 'filled');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Sensor Locations');
            axis equal;
            grid on;
        end
        function plotPotential(obj)
            % Plot the potentials at sensor locations
            scatter3(obj.Locations(:, 1), obj.Locations(:, 2), obj.Locations(:, 3), 50, obj.MeasuredV, 'filled');
            colorbar;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Sensor Potentials');
            axis equal;
            grid on;
        end
        function plotPotentialHeatmap(obj)
            % plotPotential Plots the potential as a heatmap on the sphere

            % Extract sensor locations and potential
            x = obj.Locations(:, 1);
            y = obj.Locations(:, 2);
            z = obj.Locations(:, 3);
            potential = obj.MeasuredV;


            % Create a triangulation object
            tri = delaunay(x, y, z);

            % Plot the potential as a heatmap
            figure;
            trisurf(tri, x, y, z, potential, 'EdgeColor', 'none');
            colormap(jet);
            colorbar;
            title('Potential Heatmap on Sphere');
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            grid on;
        end

        function plotEField(obj)
            % Plot the electric field vectors at sensor locations
            quiver3(obj.Locations(:, 1), obj.Locations(:, 2), obj.Locations(:, 3), ...
                    obj.MeasuredE(:, 1), obj.MeasuredE(:, 2), obj.MeasuredE(:, 3), 'AutoScale', 'on');
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Electric Field Vectors');
            axis equal;
            grid on;
        end
    end
    methods (Static)
        function locations = generateGrid(radius, thetaRange, phiRange, numTheta, numPhi, rRange, numR)
            % GENERATEGRID Generate a grid of sensor positions on a sphere.
            %
            % Inputs:
            %   radius     - Scalar radius of the sphere (used if rRange is not specified).
            %   thetaRange - 2-element vector specifying the range of polar angles [theta_min, theta_max].
            %   phiRange   - 2-element vector specifying the range of azimuthal angles [phi_min, phi_max].
            %   numTheta   - Number of points in the theta dimension.
            %   numPhi     - Number of points in the phi dimension.
            %   rRange     - (Optional) 2-element vector specifying the range of radii [r_min, r_max].
            %   numR       - (Optional) Number of points in the radial dimension.
            %
            % Output:
            %   locations  - Nx3 matrix of [x, y, z] coordinates.

            % Check if rRange and numR are specified
            if nargin < 6 || isempty(rRange)
                % If rRange is not provided, use the scalar radius
                rRange = [radius, radius]; % Single radius
                numR = 1; % Only one radius
            elseif nargin < 7
                error('If rRange is specified, numR must also be provided.');
            end

            % Generate the ranges for r, theta, and phi
            r = linspace(rRange(1), rRange(2), numR);
            theta = linspace(thetaRange(1), thetaRange(2), numTheta);
            phi = linspace(phiRange(1), phiRange(2), numPhi);

            % Create a meshgrid for the spherical coordinates
            [R, Theta, Phi] = meshgrid(r, theta, phi);

            % Convert spherical coordinates to Cartesian coordinates
            x = R .* sin(Theta) .* cos(Phi);
            y = R .* sin(Theta) .* sin(Phi);
            z = R .* cos(Theta);

            % Combine into Nx3 array of [x, y, z] locations
            locations = [x(:), y(:), z(:)];
        end

        function locations = generateRandom(radius, thetaRange, phiRange, numSensors, seed)
            % Generate random sensor positions on a sphere
            rng(seed);

            theta = randInterval(numSensors,1,thetaRange); % Random theta in thetaRange
            phi = randInterval(numSensors, 1, phiRange); % Random phi in phiRange
            r = ones(numSensors,1) * radius;

            [x,y,z] = toCartesian([r, theta, phi]); 
            locations = [x,y,z]; % Update Locations
        end
    end
end
