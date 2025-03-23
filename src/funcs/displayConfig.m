function displayConfig(config)
% Display all relevant configuration information
fprintf('--- Sensor Arrangement Configuration ---\n');
fprintf('Headshape Radius: %.2f\n', config.headshape.radius);
fprintf('Sensor Arrangement: %s\n', config.sensorArray.arrangement);
fprintf('Theta Range: [%.2f, %.2f]\n', config.sensorArray.thetaRange(1), config.sensorArray.thetaRange(2));
fprintf('Phi Range: [%.2f, %.2f]\n', config.sensorArray.phiRange(1), config.sensorArray.phiRange(2));
fprintf('Number of Theta Divisions: %d\n', config.sensorArray.numTheta);
fprintf('Number of Phi Divisions: %d\n', config.sensorArray.numPhi);
fprintf('Number of Sensors: %d\n', config.sensorArray.numSensors);
if isfield(config, 'seed')
    fprintf('Random Seed: %d\n', config.seed);
end
fprintf('-----------------------------------------\n');
end