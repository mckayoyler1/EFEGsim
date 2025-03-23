classdef TestSensorArrangement < matlab.unittest.TestCase
    methods (Test)
        function testGenerateGrid(testCase)
            % Test grid generation
            radius = 1;
            thetaRange = [0, pi];
            phiRange = [0, 2*pi];
            numTheta = 10;
            numPhi = 10;

            locations = SensorArrangement.generateGrid(radius, thetaRange, phiRange, numTheta, numPhi);
            testCase.verifySize(locations, [numTheta * numPhi, 3], ...
                'Generated grid does not have the expected size.');
        end

        function testGenerateRandom(testCase)
            % Test random location generation
            radius = 1;
            thetaRange = [0, pi];
            phiRange = [0, 2*pi];
            numSensors = 50;
            seed = 42;

            locations = SensorArrangement.generateRandom(radius, thetaRange, phiRange, numSensors, seed);
            testCase.verifySize(locations, [numSensors, 3], ...
                'Generated random locations do not have the expected size.');
        end

        function testFieldGradient(testCase)
            % Test gradient computation
            locations = [0, 0, 1; 1, 0, 0; 0, 1, 0];
            potentials = [1; 0.5; -0.2];

            fieldCalc = FieldCalculations(locations, potentials);
            [Fx, Fy, Fz] = fieldCalc.computeGradient();

            % Verify the sizes of the gradient outputs
            testCase.verifySize(Fx, [3, 1], 'Gradient Fx size is incorrect.');
            testCase.verifySize(Fy, [3, 1], 'Gradient Fy size is incorrect.');
            testCase.verifySize(Fz, [3, 1], 'Gradient Fz size is incorrect.');
        end
    end
end