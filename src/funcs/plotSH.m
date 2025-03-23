function plotSH(points, Y_vals, plotType)
% plotSH Plots precomputed spherical harmonic values on the unit sphere.
%
%   plotSH(Y_vals, points, plotType) plots the spherical harmonic values
%   provided in Y_vals at the directions specified by points.
%
%   Inputs:
%       Y_vals   - A Kx1 vector of spherical harmonic values (one harmonic)
%       points   - A Kx3 matrix where each row is [r, phi, theta] in radians.
%                  * phi: azimuth (0 to 2*pi)
%                  * theta: polar (inclination) angle (0 to pi)
%       plotType - (Optional) A string specifying what to plot:
%                  'real'   - plot the real part (default)
%                  'abs'  - plot the modulus squared, |Y|^2
%
%   Example:
%       % Suppose you have a grid of directions:
%       nPhi = 50; nTheta = 50;
%       phiVals = linspace(0, 2*pi, nPhi);
%       thetaVals = linspace(0, pi, nTheta);
%       [Phi, Theta] = meshgrid(phiVals, thetaVals);
%       pts = [Phi(:), Theta(:)];
%
%       % And you computed the spherical harmonic matrix (using your getSHMatrix):
%       Y_matrix = getSHMatrix(pts, 2); % For order N=2, columns correspond to
%                                       % (0,0), (1,-1), (1,0), (1,1), ...
%       % Choose, for example, the 4th column (which might correspond to Y_1^1)
%       Y_vals = Y_matrix(:,4);
%
%       % Now plot the real part:
%       plotSHValues(Y_vals, pts, 'real');
%
%       % Or plot the modulus squared:
%       plotSHValues(Y_vals, pts, 'abs');
%

    if nargin < 3 || isempty(plotType)
        plotType = 'real';
    end

    switch lower(plotType)
        case 'real'
            R = real(Y_vals);
        case 'abs'
            R = abs(Y_vals).^2;
        otherwise
            error('plotType must be ''real'' or ''abs''.');
    end

    % Extract phi and theta from points.
    theta = points(:,2);
    phi = points(:,3);

    % Convert spherical coordinates (phi, theta) on the unit sphere to Cartesian.
    % Here, theta is the polar (inclination) angle measured from the positive z-axis.
    [x,y,z] = toCartesian([R, theta, phi]);

    % Choose what to display.

    % Create a scatter3 plot on the unit sphere.
    figure;
    scatter3(x, y, z, 36, R, 'filled');
    axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Spherical Harmonic (%s)', plotType));
    colorbar;
    grid on;
end
