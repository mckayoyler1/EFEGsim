% EFEG Spherical Harmonic Transform
% Initialization
clc; clearvars; close all;
config = getConfig(); 
expansionOrder = 6;
basis = 'complex';
% Dipoles
dipoleDist = Dipole.generateDipoles(config);

% Dipole Potentials
syms x y z
symbV = computeSymbolicPotential(dipoleDist,x,y,z);
symbE = computeEFromV(symbV, x, y, z);
% Data
% Make Grid Points
numTheta = 50;
numPhi = 120;
gridPos = SensorArray.generateGrid(config.headshape.radius, [pi/(2*numTheta),pi/2],[2*pi/numPhi, 2*pi], numTheta, numPhi);
randPos = SensorArray.generateRandom(config.headshape.radius, [pi/(2*numTheta),pi/2],[2*pi/numPhi, 2*pi], numTheta*numPhi, config.seed);

% Compute Spherical Harmonic Expansion
dataRand = computeSHExpansion(randPos, evaluateSymV(randPos, symbV),expansionOrder, basis);
dataGrid = computeSHExpansion(gridPos, evaluateSymV(gridPos, symbV),expansionOrder, basis);
% Compute E Fields
dataRand = getEfromSHcoeffs(dataRand);
dataRand.E.trueE = evaluateSymE(symbE, dataRand.points);

dataGrid = getEfromSHcoeffs(dataGrid);
dataGrid.E.trueE = evaluateSymE(symbE, dataGrid.points);

% TEST DATA
%%%%% TEST DATA %%%%%%
Y_N_Grid = getSH2L(9, dataGrid.sphPoints(:,2:3), basis);
Y_N_Rand = getSH2L(9, dataRand.sphPoints(:,2:3), basis);
dataRandTest = computeSHExpansion(randPos, real(getSHEntry(Y_N_Rand, 3,2)), expansionOrder, basis); % V is spherical harmonic Y_(3,2)
dataGridTest = computeSHExpansion(gridPos, real(getSHEntry(Y_N_Grid, 3,2)), expansionOrder, basis);

dataRandTest = getEfromSHcoeffs(dataRandTest);
dataRandTest.E.trueE = evaluateSymE(symbE, dataRandTest.points);

dataGridTest = getEfromSHcoeffs(dataGridTest);
dataGridTest.E.trueE = evaluateSymE(symbE, dataGridTest.points);


% Visualize
data = dataRand;
figure
%%% Plot Potentials
subplot(2,3,1), ft_plot_topo3d(data.points, data.V.trueV, 'neighbourdist', inf); title('true V');
colorbar;
subplot(2,3,2), ft_plot_topo3d(data.points, data.V.reconstructedV, 'neighbourdist', inf); title('Reconstructed V');
colorbar;
subplot(2,3,3), ft_plot_topo3d(data.points, abs(data.V.reconstructedV - data.V.trueV), 'neighbourdist', inf); title('Error in V');
colorbar;

figure
%%% Plot E Fields
plotSphericalVectorField(data.sphPoints, data.E.sphReal); title('Reconstructed E Field');
figure
axis equal;
plotCartesianVectorField(data.points, data.E.trueE); title('True E Field');
figure
axis equal;
plotCartesianVectorField(data.points,data.E.trueE - data.E.cartReal); title('Error in E field'); title('Error in E');
axis equal;
%data.viz.EFieldPlots = EFieldPlots;
%vizSHExpansion(dataRand);
%vizSHExpansion(dataGridTest);


%% Check V With Package
Nord = expansionOrder;
interpGrid = dataGrid.points;
interpSPHGrid = toSphere(interpGrid);
SPHang = [interpSPHGrid(:,3), interpSPHGrid(:,2)];
SPHang(:,2) = pi/2 - SPHang(:,2);
F_reg = dataGrid.V.trueV;


voronoi_Weights = getVoronoiWeights(SPHang);

Fnm_reg = leastSquaresSHT(Nord, F_reg, SPHang, 'complex', voronoi_Weights);
Finterp = inverseSHT(Fnm_reg, SPHang, 'complex');
Vtest = real(Finterp);
figure;
ft_plot_topo3d(dataGrid.points, Vtest, 'neighbourdist', inf);
colorbar;
figure;
ft_plot_topo3d(dataGrid.points, dataGrid.V.trueV - Vtest, 'neighbourdist', inf); title('error in V from Package')
colorbar;
figure;
ft_plot_topo3d(dataGrid.points, dataGrid.V.trueV, 'neighbourdist', inf); title('True V')
colorbar;

