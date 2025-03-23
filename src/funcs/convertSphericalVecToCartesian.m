function cartVectors = convertSphericalVecToCartesian(sphPoints, sphVectors)
% convertSphericalToCartesian Converts spherical points and vectors to Cartesian.
%
%   [cartPoints, cartVectors] = convertSphericalToCartesian(sphPoints, sphVectors)
%
%   INPUTS:
%       sphPoints  - An Mx3 matrix where each row is [r, theta, phi] (angles in radians)
%       sphVectors - An Mx3 matrix where each row is [V_r, V_theta, V_phi]
%
%   OUTPUTS:
%       cartPoints  - An Mx3 matrix of Cartesian coordinates [x, y, z]
%       cartVectors - An Mx3 matrix of Cartesian vector components [V_x, V_y, V_z]
%
%   The conversion formulas are:
%
%       x = r*sin(theta)*cos(phi)
%       y = r*sin(theta)*sin(phi)
%       z = r*cos(theta)
%
%       V_x = V_r*sin(theta)*cos(phi) + V_theta*cos(theta)*cos(phi) - V_phi*sin(phi)
%       V_y = V_r*sin(theta)*sin(phi) + V_theta*cos(theta)*sin(phi) + V_phi*cos(phi)
%       V_z = V_r*cos(theta) - V_theta*sin(theta)
%
%   Example:
%       % Define some spherical points and vectors:
%       sphPoints = [1, pi/4, pi/4; 1, pi/3, pi/2];  % each row is [r, theta, phi]
%       sphVectors = [1, 0, 0; 0, 1, 0];              % example vector components
%
%       [cartPts, cartVecs] = convertSphericalToCartesian(sphPoints, sphVectors);

    % Extract spherical coordinates.
    r     = sphPoints(:,1);
    theta = sphPoints(:,2);
    phi   = sphPoints(:,3);
    
    % Convert points.
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
    
    
    % Extract spherical vector components.
    Vr     = sphVectors(:,1);
    Vtheta = sphVectors(:,2);
    Vphi   = sphVectors(:,3);
    
    % Convert vector components.
    Vx = Vr .* sin(theta) .* cos(phi) + Vtheta .* cos(theta) .* cos(phi) - Vphi .* sin(phi);
    Vy = Vr .* sin(theta) .* sin(phi) + Vtheta .* cos(theta) .* sin(phi) + Vphi .* cos(phi);
    Vz = Vr .* cos(theta) - Vtheta .* sin(theta);
    
    cartVectors = [Vx, Vy, Vz];
end