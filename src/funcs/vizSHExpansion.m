function data = vizSHExpansion(data)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


%%% Plot Potentials
Vplots = figure('Name', 'Potential Plots','Color','w','Position',[100 100 900 700]);
subplot(3,1,1), ft_plot_topo3d(data.points, data.V.trueV, 'neighbourdist', inf); title('true V');
colorbar;
subplot(3,1,2), ft_plot_topo3d(data.points, data.V.reconstructedV, 'neighbourdist', inf); title('Reconstructed V');
colorbar;
subplot(3,1,3), ft_plot_topo3d(data.points, abs(data.V.reconstructedV - data.V.trueV), 'neighbourdist', inf); title('Error in V');
colorbar;
zoom()
data.viz.Vplots = Vplots;

%%% Plot E Fields
EFieldPlots = figure('Name', 'E Vector Field Plots','Color','w','Position',[100 100 900 700]);
subplot(1,3,1), plotSphericalVectorField(data.sphPoints, data.E.sphReal); title('Reconstructed E Field');
axis equal;
subplot(1,3,2), plotCartesianVectorField(data.points, data.E.trueE); title('True E Field');
axis equal;
subplot(1,3,3), plotCartesianVectorField(data.points,data.E.trueE - data.E.cartReal); title('Error in E field'); title('Error in E');
axis equal;
data.viz.EFieldPlots = EFieldPlots;

%EFieldErrorPlots = figure('Name', 'E Field Error Visualizations');

end