function ftData = toFieldTripData(sensorArray)
    % Convert SensorArray potentials to FieldTrip-compatible data struct
    ftData.avg = sensorArray.MeasuredV;   % Sensor potential values (Nx1 vector)
    ftData.time = 0;                       % Time (e.g., 0 for static field)
    ftData.dimord = 'chan_time';           % Specify data organization
end