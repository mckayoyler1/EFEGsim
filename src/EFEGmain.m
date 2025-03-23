%% Initialization
clc; clearvars; close all;
config = getConfig(); 

% Generate dipoles
dipoleList = Dipole.generateDipoles(config);

testDipole = Dipole([0,4,4],[1,2,1]);
dipoles = dipoleList;

syms x y z
symbolicV = computeSymbolicPotential(dipoles,x,y,z);
symbolicE = computeEFromV(symbolicV, x, y, z);


grid40t120pcfg = struct();
grid40t120pcfg.radius = config.headshape.radius;
grid40t120pcfg.Ntheta = 40;
grid40t120pcfg.Nphi   = 120;

interpGrid = genSphereGrid(grid40t120pcfg);


% Generate Grid Sensor Arrangement
sensorArrayGrid = SensorArray(config);
sensorArrayGrid.measurePotential(dipoles);






































% RBF Interp
Vx = rbf_interpolation_spherical(sensorArrayGrid.Locations, sensorArrayGrid.MeasuredV, interpGrid, 'gaussian', .7);
Vfunc = matlabFunction(symbolicV, 'Vars', [x, y, z]);
Vtrue = Vfunc(interpGrid(:, 1), interpGrid(:, 2), interpGrid(:, 3));
Vdiff = abs((Vx - Vtrue));



figure('Name','rbf')
scatter3(interpGrid(:, 1), interpGrid(:, 2), interpGrid(:, 3), 50, Vdiff, 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('RBF Potentials');
axis equal;
grid on;

figure('Name','Slices')
hold on;
p1 = plotSlice(interpGrid,Vtrue);
p1.Color = "red";
hold on
p2 = plotSlice(interpGrid, Vx);
p2.Color = "blue";


%%
% config.sensorArray.arrangement = 'random';
% sensorArrayRand = SensorArray(config);
% sensorArrayRand.measurePotential(dipoles);
% 
% plotInterpV_ft(sensorArrayGrid)
% plotDipoles(dipoles)
% %sensorArrayGrid.plotSensors()
% 
% sensorArrayGrid.plotPotential()
% plotSymbolicEField(symbolicE, sensorArrayGrid.Locations);
%sensorArrayGrid.plotEField()


% VxMesh = unflattenMeshgrid(interpGrid, Vx, numTheta, numPhi);
% % VxMesh (row decides theta, column decides phi, [1,2,3] picks
% % [r,theta,phi] component so V(r,theta,phi) = VxMesh(:,:,1:3)




Etrue = evaluateSymE(symbolicE, interpGrid);
gradientsx = gradient_spherical(interpGrid,Vx);
Ex = -gradientsx;
cartEx = -toCartesianGradient(gradientsx, interpGrid);
Ediff = Etrue-cartEx;
figure('Name','E from Interpolated V')
plotCartesianVectorField(interpGrid, Ex)
title('E from Interpolated V')
figure('Name','True E')
plotCartesianVectorField(interpGrid, Etrue)

figure('Name', 'Difference of Etrue and E interpolated')
plotCartesianVectorField(interpGrid, Ediff)
displayConfig(config)