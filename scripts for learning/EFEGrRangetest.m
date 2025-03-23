clc; clearvars; close all;
config = getConfig(); 

% Generate dipoles
dipoleList = Dipole.generateDipoles(config);

testDipole = Dipole([0,0,8],[0,0,4]);
dipoles = dipoleList;

syms x y z
symbolicV = computeSymbolicPotential(dipoles,x,y,z);
symbolicE = computeEFromV(symbolicV, x, y, z);
numTheta = 40;
numPhi = 120;
numR = 3;
interpGrid = SensorArray.generateGrid( ...
    config.headshape.radius, ...
    [pi/(2*numTheta),pi/2],[2*pi/numPhi, 2*pi], ...
    numTheta, numPhi, config.headshape.rRange, numR);

interpGridScalp = SensorArray.generateGrid( ...
    config.headshape.radius, ...
    [pi/(2*numTheta),pi/2],[2*pi/numPhi, 2*pi], ...
    numTheta, numPhi);

% Generate Grid Sensor Arrangement
sensorArrayGrid = SensorArray(config);
sensorArrayGrid.measurePotential(dipoles);


%% RBF Interp
Vx = rbf_interpolation_spherical(sensorArrayGrid.Locations, sensorArrayGrid.MeasuredV, interpGrid, 'gaussian', .35);
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
hold on
dimscalp = size(interpGridScalp);
length = dimscalp(1);
p3 = plotSlice(interpGridScalp, Vtrue(1:length));
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
% VxMesh (row decides theta, column decides phi, [1,2,3] picks
% [r,theta,phi] component so V(r,theta,phi) = VxMesh(:,:,1:3)



Etrue = evaluateSymE(symbolicE, interpGrid);
gradientsx = gradient_spherical(interpGrid,Vx);
Ex = -gradientsx;
cartEx = -toCartesianGradient(gradientsx, interpGrid);
Ediff = Etrue-cartEx;

figure('Name','E from Interpolated V')
plot_EfieldScalp(interpGrid, Ex, Vx, 'spherical')

figure('Name','True E')
plot_EfieldScalp(interpGrid, Etrue, Vtrue, 'cartesian')

figure('Name', 'Difference of Etrue and E interpolated')
plot_EfieldScalp(interpGrid, Ediff, Vdiff)
displayConfig(config)