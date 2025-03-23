function vizSH(Y, valueType, scaleFactor)
% vizSH Visualizes a symbolic spherical harmonic as an orbital-like surface.
%
%   vizSH(Y) plots the spherical harmonic Y as an orbital by modulating
%   the radius of a spherical grid using the harmonic’s value.
%
%   vizSH(Y, valueType) lets you choose how to modulate the radius:
%         'real' - use the real part of Y
%         'abs'  - use the absolute value of Y (default)
%
%   vizSH(Y, valueType, scaleFactor) applies an additional scaling factor to
%   the radial deformation. (Default is 1.)
%
%   INPUT:
%       Y          - A symbolic expression for the spherical harmonic in terms 
%                    of theta and phi.
%       valueType  - (Optional) 'real' or 'abs' (default is 'abs').
%       scaleFactor- (Optional) A scalar factor for the deformation (default = 1).
%
%   Example:
%       syms theta phi
%       Y = symbolicYlm(2, 1);  % your symbolic spherical harmonic
%       vizSH(Y, 'real', 0.5);

    %% Set default values for optional arguments
    if nargin < 2 || isempty(valueType)
        valueType = 'abs';
    end
    if nargin < 3 || isempty(scaleFactor)
        scaleFactor = 1;
    end

    %% Ensure symbolic variables exist and convert Y to a function handle
    syms theta phi
    try
        f = matlabFunction(Y, 'Vars', [theta, phi]);
    catch ME
        error('Failed to convert Y to a function handle. Ensure Y is expressed in terms of theta and phi.');
    end

    %% Create a grid of spherical angles using ndgrid.
    nTheta = 150; % resolution for theta (polar angle)
    nPhi = 150;   % resolution for phi (azimuth)
    thetaVals = linspace(0, pi, nTheta);
    phiVals   = linspace(0, 2*pi, nPhi);
    [Theta, Phi] = ndgrid(thetaVals, phiVals);
    
    %% Evaluate the spherical harmonic on the angular grid.
    Y_vals = f(Theta, Phi);
    
    %% Determine the modulation to apply to the radius.
    switch lower(valueType)
        case 'real'
            modulation = real(Y_vals);
        case 'abs'
            modulation = abs(Y_vals);
        otherwise
            error('Unknown valueType. Choose either ''real'' or ''abs''.');
    end

    %% Compute the modulated radius.
    % Start with the unit sphere (radius = 1) and add a deformation.
    R = 1+ scaleFactor * modulation;

    %% Convert the modulated spherical coordinates back to Cartesian coordinates.
    % Here theta is the polar angle and phi is the azimuth.
    X = R .* sin(Theta) .* cos(Phi);
    Y_cart = R .* sin(Theta) .* sin(Phi);
    Z = R .* cos(Theta);

    %% Plot the orbital-like surface.
    figure;
    % If you prefer a uniform color, remove the fourth argument from surf.
    h = surf(X, Y_cart, Z, 'EdgeColor', 'none'); 
    colormap jet;
    axis equal;
    shading interp;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Orbital Visualization (%s modulation)', valueType));
    
    % Add lighting for a 3D effect.
    camlight headlight;
    lighting phong;
    material shiny;
    view(40,30);
end
