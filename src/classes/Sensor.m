classdef Sensor
    properties
        Location  % Sensor location as a 1x3 vector [x, y, z]
        Potential % Measured potential at the sensor (scalar)
        EField
        Label
    end

    methods
        function obj = Sensor(location, potential)
            if nargin < 1
                location = [0, 0, 0]; % Default location
            end
            % Constructor to initialize the sensor
            if nargin < 2
                potential = 0; % Default potential if not provided
            end
            obj.Location = location;
            obj.Potential = potential;
            obj.EField = [];
        end

        function obj = measurePotential(obj, dipoles)
            % Measure the potential at the sensor due to multiple dipoles
            % Input:
            %   dipoles - Array of Dipole objects
            %   k       - Scalar constant (e.g., proportionality constant)

            numDipoles = numel(dipoles);
            totalPotential = 0;

            for i = 1:numDipoles
                totalPotential = totalPotential + dipoles(i).getPotential(obj.Location);
            end

            obj.Potential = totalPotential;
            obj.Label = 'sensorX';
        end
    end

    methods (Static)
        function sensors = generateSensors(sensorArray)
            % Generate an array of sensors from a set of locations
            % Input:
            %   sensorArrangement object
            % Output:
            %   sensors   - Array of Sensor objects
            locations = sensorArray.Locations;
            numSensors = sensorArray.NumSensors;
            sensors(numSensors, 1) = Sensor(); % Preallocate for efficiency

            for i = 1:numSensors
                sensors(i) = Sensor(locations(i, :));
                sensors(i).Label = sprintf('Sensor%d', i);
            end
        end
    end
end