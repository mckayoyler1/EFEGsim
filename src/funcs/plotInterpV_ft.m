function plotInterpV_ft(sensorArray)

%ft_topoplotER(cfg, data);

%plotSensors(sensorArray)
ft_plot_topo3d(sensorArray.Locations, sensorArray.MeasuredV, 'neighbourdist', inf);
view(3);
end