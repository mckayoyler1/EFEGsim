function plotV(points, V)

%ft_topoplotER(cfg, data);

%plotSensors(sensorArray)
ft_plot_topo3d(points,V, 'neighbourdist', inf);
view(3);
end