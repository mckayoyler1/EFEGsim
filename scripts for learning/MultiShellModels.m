%%%%%%% YURY PETROV E FIELD AND POTENTIAL %%%%%%%%%%
clc; clearvars; close all;
% Example model setup with 4 shells (brain, CSF, skull, scalp)
% model_params.epsilon1 = 0.3;  % brain conductivity (S/m)
% 
% % Define shell parameters
% model_params.shells(1).inner_radius = 0;
% model_params.shells(1).outer_radius = 0.79;  % 7.9 cm
% model_params.shells(1).epsilon = 0.3;  % radial conductivity
% model_params.shells(1).eta = 0.3;      % tangential conductivity
% 
% model_params.shells(2).inner_radius = 0.79;
% model_params.shells(2).outer_radius = 0.8;  % 8.0 cm
% model_params.shells(2).epsilon = 1.5;  % CSF radial conductivity
% model_params.shells(2).eta = 1.5;      % CSF tangential conductivity
% 
% model_params.shells(3).inner_radius = 0.8;
% model_params.shells(3).outer_radius = 0.85;  % 8.5 cm
% model_params.shells(3).epsilon = 0.01;  % Skull radial conductivity
% model_params.shells(3).eta = 0.023;     % Skull tangential conductivity
% 
% model_params.shells(4).inner_radius = 0.85;
% model_params.shells(4).outer_radius = Inf;    % Scalp extends to infinity
% model_params.shells(4).epsilon = 0.3;  % Scalp radial conductivity
% model_params.shells(4).eta = 0.3;      % Scalp tangential conductivity
% 
% % Example usage
% r0 = .7;              % 7 cm source location
% alpha = 0;            % radial dipole
% q = 3e-4;             % dipole magnitude (A·m)

config = getConfig(); 
data = struct('expansionOrder', 40        , ...
              'basis'         , 'complex', ...
              'cfg'           , config   ,...
              'Dipoles'       , Dipole.generateDipoles(config), ...
              'V'             , struct(), ...
              'E'             , struct(), ...
              'grid_type'     , 'uniform');
data = getShellPoints(data);

% potential = calculate_dipole_potential(data.sphpos, r0, alpha, q, model_params);
% data.V.true = potential;
% data.E.true = calculate_dipole_electric_field(data.sphpos, r0, alpha, q, model_params);
%%
data = computeSHExpansion(data);

%%
Escale = 1e6;
plotdata(data,Escale)
% figure;
% plotSphericalVectorField(data.sphpos, data.E.true, 1e6); title('True E Field');
% axis equal;
% 
% figure('Name','True V');
% ft_plot_topo3d(data.pos,real(data.V.true), 'neighbourdist', inf);
% title('True V');
% colorbar;
% view(3);

%%

function plotdata(data, Escale)
    if nargin < 2 || isempty(Escale)
        Escale = 1;
    end
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
plotSphericalVectorField(data.sphpos, data.E.sph, Escale); title('Reconstructed E Field');
axis equal;

% figure;
% plotCartesianVectorField(data.pos, data.E.cart); title('Cartesian E Field');
% axis equal;

figure;
plotSphericalVectorField(data.sphpos, data.E.true, Escale); title('True E Field');
axis equal;

figure;
plotSphericalVectorField(data.sphpos, (data.E.sph - data.E.true), Escale); title('Scaled Difference (E_{recon} - E_{true})*1e6');
axis equal;
end