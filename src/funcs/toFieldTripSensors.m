function elec = toFieldTripSensors(sensorArray)
    % Convert SensorArray object to FieldTrip-compatible sensor struct
    elec.chanpos = sensorArray.Locations;   % Sensor positions (Nx3 matrix)
    elec.label = arrayfun(@(i) sprintf('Sensor%d', i), 1:size(sensorArray.Locations, 1), 'UniformOutput', false);
    elec.unit = 'cm'; % Or 'mm' depending on your scale
end