%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TROUBLESHOOTING SH EXPANSION %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%PLEASE 

clc; clearvars; close all;
config = getConfig(); 
syms x y z
% 
% data = struct('expansionOrder', 40        , ...
%               'basis'         , 'complex', ...
%               'cfg'           , config   ,...
%               'Dipoles'       , Dipole.generateDipoles(config), ...
%               'V'             , struct(), ...
%               'E'             , struct(), ...
%               'grid_type'     , 'modulated');
% data.V.symbV = computeSymbolicPotential(data.Dipoles, x, y, z);
% data.E.symbE = computeEFromV(data.V.symbV, x, y, z);
% data = getShellPoints(data);
% data = computeSHExpansion(data);

%veronoiVerify = verifyVoronoiWeights(data);
%%
%plotdata(data) % Plots modulated data

%%
%%%%%%%% UNIFORM RADIUS GRID %%%%%%%%%%%
data = struct('expansionOrder', 40        , ...
              'basis'         , 'complex', ...
              'cfg'           , config   ,...
              'Dipoles'       , Dipole.generateDipoles(config), ...
              'V'             , struct(), ...
              'E'             , struct(), ...
              'grid_type'     , 'uniform');
data.V.symbV = computeSymbolicPotential(data.Dipoles, x, y, z);
data.E.symbE = computeEFromV(data.V.symbV, x, y, z);
data = getShellPoints(data);
data = computeSHExpansion(data);
%%
%veronoiVerify = verifyVoronoiWeights(data);
plotdata(data)

%%


function plotdata(data)
combined_min = min(min(data.V.true), min(real(data.V.complexSHT)));
combined_max = max(max(data.V.true), max(real(data.V.complexSHT)));
cscale = real([combined_min combined_max]);
figure('Name','Reconstructed V');
ft_plot_topo3d(data.pos,real(data.V.complexSHT), 'neighbourdist', inf);
colorbar;
clim(cscale);
view(3);
title('Reconstructed V');

figure('Name','True V');
ft_plot_topo3d(data.pos,real(data.V.true), 'neighbourdist', inf);
title('True V');
colorbar;
clim(cscale);
view(3);

figure('Name','Error in V');
ft_plot_topo3d(data.pos,real(data.V.true-data.V.complexSHT), 'neighbourdist', inf);
title('Error in V');
colorbar;
view(3);

figure
scatter3(data.pos(:,1),data.pos(:,2),data.pos(:,3),36,real(data.V.complexSHT), 'filled');
colorbar;
clim(cscale);
axis equal;
title('Scatter plot of Reconstructed V')
view(3);

figure
scatter3(data.pos(:,1),data.pos(:,2),data.pos(:,3),36,real(data.V.true), 'filled');
colorbar;
clim(cscale);
axis equal;
title('Scatter plot of true V')
view(3);

figure
scatter3(data.pos(:,1),data.pos(:,2),data.pos(:,3),36,real(data.V.true - data.V.complexSHT), 'filled');
colorbar;
axis equal;
title('Scatter plot of V error')
view(3);
% 
% plotSH(data.sphpos,real(data.V.true-data.V.complexSHT),'abs')
% axis equal;

figure;
plotSphericalVectorField(data.sphpos, data.E.sph); title('Reconstructed E Field');
axis equal;

% figure;
% plotCartesianVectorField(data.pos, data.E.cart); title('Cartesian E Field');
% axis equal;

figure;
plotCartesianVectorField(data.pos, data.E.true); title('True E Field');
axis equal;
end