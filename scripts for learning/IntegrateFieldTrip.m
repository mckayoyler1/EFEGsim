clc; clearvars;
config = getConfig(); 

% Generate dipoles
dipoleList = Dipole.generateDipoles(config);

testDipole = Dipole([0,0,.08],[0,0,4]);
dipoles = dipoleList;

sensorArray = SensorArray(config);
sensorArray.measurePotential(dipoles);

data = struct();
data.avg = sensorArray.MeasuredV;

labels = cell(1, sensorArray.NumSensors); % Preallocate
for i = 1:sensorArray.NumSensors
    labels{i} = sprintf('Sensor%d', i);
end
data.label = labels;
data.time = 0;
data.mord = 'chan_time';

% Convert to a FieldTrip layout structure
layout = struct();
layout.pos = sensorArray.Locations(:, 1:2); % Use x, y coordinates only
layout.label = labels; % Sensor labels
layout.width = 0.03 * ones(size(labels)); % Dummy width for each sensor
layout.height = 0.03 * ones(size(labels)); % Dummy height for each sensor

% Configure and plot
cfg = [];
cfg.layout = layout; % Use custom layout
cfg.parameter = 'avg'; % Specify the parameter to plot
%ft_topoplotER(cfg, data);

%plotSensors(sensorArray)
ft_plot_topo3d(sensorArray.Locations, sensorArray.MeasuredV, 'neighbourdist', inf);