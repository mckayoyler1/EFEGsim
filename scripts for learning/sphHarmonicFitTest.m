clc; clearvars; close all;
config = getConfig(); 


% % Generate dipoles
dipoleList = Dipole.generateDipoles(config);
% 
% testDipole = Dipole([0,0,.8],[0,0,5]);
dipoles = dipoleList;
% 
syms x y z
symbolicV = computeSymbolicPotential(dipoles,x,y,z);
symbolicE = computeEFromV(symbolicV, x, y, z);
numTheta = 40;
numPhi = 120;
interpGrid = SensorArray.generateGrid( ...
    config.headshape.radius, ...
    [pi/(2*numTheta),pi],[2*pi/numPhi, 2*pi], ...
    numTheta, numPhi);
Vfunc = matlabFunction(symbolicV, 'Vars', [x, y, z]);
Vtrue = Vfunc(interpGrid(:, 1), interpGrid(:, 2), interpGrid(:, 3));


% % Generate Grid Sensor Arrangement
% sensorArrayGrid = SensorArray(config);
% % sensorArrayGrid.measurePotential(dipoles);
% sensorArrayGrid.MeasuredV = Vfunc(sensorArrayGrid.Locations(:, 1), sensorArrayGrid.Locations(:, 2), sensorArrayGrid.Locations(:, 3));

interpSPHGrid = toSphere(interpGrid);
SPHang = [interpSPHGrid(:,3), interpSPHGrid(:,2)];
SPHang(:,2) = pi/2 - SPHang(:,2);


% testY = getSH(4, SPHang, 'complex');
% test2 = getSH2L(10, interpSPHGrid(:,2:3), 'complex');
% test3 = spharm(interpSPHGrid(:,2), interpSPHGrid(:,3), 3,2);

%%
 voronoiWeights = getVoronoiWeights(interpSPHGrid(:,2:3));
 
 %%

[x,y,z] = toCartesian([abs(real(test3)),interpSPHGrid(:,2), interpSPHGrid(:,3)]);
%scatter3(x, y, z, 50,abs(real(test2(:,8))) , 'filled')
% 
% plotSH(interpSPHGrid,test2(:,3),  'abs')
scatter3(x,y,z,50, abs(real(test3)))

%%
result = computeSHExpansion(interpGrid, Vtrue, 30, 'complex');
Etest = getEfromSHcoeffs(result);



V_SPH = result.V.reconstructedV;
figure
subplot(121), ft_plot_topo3d(interpGrid, Vtrue, 'neighbourdist', inf); title('true V')

subplot(122), ft_plot_topo3d(interpGrid, V_SPH, 'neighbourdist', inf); title('Reconstructed V')
colorbar;

figure('Name',  'Difference')
ft_plot_topo3d(interpGrid, abs((Vtrue-V_SPH)./Vtrue), 'neighbourdist', inf); title('Difference')
colorbar;

Nord = 15;
Y11_example_fit = leastSquaresSHT(Nord, testY(:,230), SPHang, 'complex', test_weights);
% plot coefficients
figure
plot(0:(Nord+1)^ 2-1, abs(Y11_example_fit))
hold on
title('SHD spectral coefficients for Y_{11}'), xlabel('q = n^2+n+m'), ylabel('Magnitude')



%%
clc; clearvars; close all;
config = getConfig(); 

dipoleList = Dipole.generateDipoles(config);
testDipole = Dipole([0,0,0.08],[0,0,3]);
dipoles = dipoleList;

syms x y z
symbV = computeSymbolicPotential(dipoles,x,y,z);
symbE = computeEFromV(symbV, x, y, z);
Vfunc = matlabFunction(symbV, 'Vars', [x, y, z]);
numTheta = 30;
numPhi = 200;
interpGrid = SensorArray.generateGrid( ...
    config.headshape.radius, ...
    [pi/(2*numTheta),pi/2],[2*pi/numPhi, 2*pi], ...
    numTheta, numPhi);
interpSPHGrid = toSphere(interpGrid);

syms theta phi real
data = struct();
Vtrue = Vfunc(interpGrid(:, 1), interpGrid(:, 2), interpGrid(:, 3));
Y_N = getSH2L(9, interpSPHGrid(:,2:3), 'complex');
data.V = struct();
data.V.trueV = real(getSHEntry(Y_N, 3,2));% getSHEntry(Y_N, 1,1) + getSHEntry(Y_N, 2,1) + getSHEntry(Y_N, 7,3));


data = computeSHExpansion(interpGrid, data.V.trueV, 9, 'complex');

% figure
% subplot(2,2,1), ft_plot_topo3d(interpGrid, data.V.trueV, 'neighbourdist', inf); title('true V')
% colorbar;
% subplot(2,2,2), ft_plot_topo3d(interpGrid, data.V.reconstructedV, 'neighbourdist', inf); title('Reconstructed V')
% colorbar;
% subplot(2,2,3), ft_plot_topo3d(interpGrid, abs(data.V.reconstructedV - data.V.trueV), 'neighbourdist', inf); title('Error in V')
% colorbar;
% figure
% ft_plot_topo3d(interpGrid, abs((data.V.myreconV-data.V.trueV)./data.V.trueV), 'neighbourdist', inf); title('Absolute Error in my recon V')
% colorbar; 

data = getEfromSHcoeffs(data);
data.E.cartesian = convertSphericalVecToCartesian(data.sphPoints, data.E.sphReal);
data.E.trueE = evaluateSymE(symbE, interpGrid);
TEST_E = -sum(data.E.test_check,3);
Test_E_cart = convertSphericalVecToCartesian(data.sphPoints, TEST_E);
Test_E_cart = real(Test_E_cart);
data.E.trueE = Test_E_cart;

vizSHExpansion(data);
% 
% figure
% subplot(1,3,1), ft_plot_topo3d(interpGrid, data.V.trueV, 'neighbourdist', inf); title('true V')
% colorbar;
% subplot(1,3,2), ft_plot_topo3d(interpGrid, data.V.reconstructedV, 'neighbourdist', inf); title('Reconstructed V')
% colorbar;
% subplot(1,3,3), ft_plot_topo3d(interpGrid, abs((data.V.reconstructedV-data.V.trueV)./data.V.trueV), 'neighbourdist', inf); title('Absolute Error in V')
% colorbar; 
% vectorCheck = Y_N*data.E.scaledCoeffs;
% data.E.fromGrad = -gradient_spherical(interpGrid,data.V.trueV);
% 
% figure;
% subplot(1,3,1), plotSphericalVectorField(data.sphPoints, data.E.sphReal); title('Reconstructed E field')
% % plotSphericalVectorField(data.sphPoints, data.E.fromGrad); title('E from Numerical Grad(V)')
% subplot(1,3,2), plotCartesianVectorField(data.points, data.E.trueE); title('True E field')
% plotCartesianVectorField(data.points, Test_E_cart); title('Test E field')
% plotCartesianVectorField(data.points,data.E.trueE - data.E.cartesian); title('Error in E field')

%%
config = getConfig(); 



numTheta = 10;
numPhi = 200;
interpGrid = SensorArray.generateGrid( ...
    config.headshape.radius, ...
    [pi/(2*numTheta),pi],[2*pi/numPhi, 2*pi], ...
    numTheta, numPhi);
interpSPHGrid = toSphere(interpGrid);

r = interpSPHGrid(:,1);
theta = interpSPHGrid(:,2);
phi = interpSPHGrid(:,3);


l = 1;
Plm = legendre(l,cos(theta)');

m = 0;
if l ~= 0
    Plm = Plm(m+1,:)';
end

a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
Exp = exp(1i*m*phi);
Ylm = C .* Plm .*Exp;

[x,y,z] = toCartesian([abs(real(Ylm)),theta, phi]);

F = scatteredInterpolant(x, y, z);
% Define a coarser grid to keep memory usage down:
[xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
zq = F(xq, yq);
surf(xq, yq, zq)

scatter3(x,y,z)
% surf(Xm,Ym,Zm)
title('$Y_3^2$ spherical harmonic','interpreter','latex')

% [Xm,Ym,Zm] = sph2cart(phi, pi/2-theta, abs(real(Ylm)));
% surf(Xm,Ym,Zm)
% title('$Y_3^2$ spherical harmonic','interpreter','latex')

%%
%%

dx = pi/60;
col = 0:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,col);
l = 4;
Plm = legendre(l,cos(theta));
m = 2;
if l ~= 0
    Plm = reshape(Plm(m+1,:,:),size(phi));
end
a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
Ylm = C .*Plm .*exp(1i*m*phi);
[Xm,Ym,Zm] = sph2cart(phi, pi/2-theta, abs(real(Ylm)));
surf(Xm,Ym,Zm)
title('$Y_3^2$ spherical harmonic','interpreter','latex')