classdef Dipole
    properties
        Position  % Position of the dipole as a 1x3 vector [x, y, z]
        Moment    % Dipole moment as a 1x3 vector
    end
    
    methods
        function obj = Dipole(position, moment)
            % Constructor to initialize the dipole
            obj.Position = position;
            obj.Moment = moment;
        end
        

        function potential = V(obj, observationPoints)
            % Method to compute the potential at a given set of points
            etas = observationPoints - obj.Position;
            etaNorm = vecnorm(etas,2,2);
            
            if etaNorm == 0
                error('Observation point cannot coincide with the dipole position.');
            end
            potential = Constants.k * sum(etas .* obj.Moment, 2) ./ (etaNorm.^3);
            
        end
    end

    methods (Static)
        function moments = generateMoments(positions, numdips, momentType, momentMagnitude)
            % Generate random or radial dipole moments. 
            % uses numDips, momentType, and Q from dipoleParams
            
            if momentType == "radial"
                x = positions(:, 1);
                y = positions(:, 2);
                z = positions(:, 3);
                norms = sqrt(x.^2 + y.^2 + z.^2);
                Q = momentMagnitude;
                moments = Q * [x./norms, y./norms, z./norms];
            elseif  momentType == "random"             
             moments = randn(numdips, 3); % Nx3 matrix of random moments
            else 
                disp("ERROR: invalid dipole moment type")
            end
        end

        function dipoles = generateDipoles(config)
            % Generate random dipole positions within specified ranges

            % Extract ranges from config
            rRange = config.dipoleParams.rRange;
            thetaRange = config.dipoleParams.thetaRange;
            phiRange = config.dipoleParams.phiRange;
            
            % Set random seed
            seed = config.seed;
            rng(seed);
            % Generate random spherical coordinates
            numdips = config.dipoleParams.numDips;
            r = randInterval(numdips, 1, rRange);
            theta = randInterval(numdips, 1, thetaRange);
            phi = randInterval(numdips, 1, phiRange);

            % Convert spherical coordinates to Cartesian
            spherical = [r, theta, phi];
            [X,Y,Z] = toCartesian(spherical);
            positions = [X,Y,Z];
            % Generate random or specific dipole moments
            momentType = config.dipoleParams.momentType;
            momentMagnitude = config.dipoleParams.momentMagnitude;
            moments = Dipole.generateMoments(positions, numdips, momentType, momentMagnitude); % A function to generate moments

            % Create an array of Dipole objects
            dipoles(1:numdips, 1) = Dipole(zeros(1,3),zeros(1,3)); % Preallocate for efficiency
            for i = 1:numdips
                dipoles(i) = Dipole(positions(i, :), moments(i, :));
            end
        end
    end
end
