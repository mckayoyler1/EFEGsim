function spherical = toSphere(cartesian)
% CARTESIANTOSPHERICAL Converts a Cartesian vector to spherical coordinates.
% Input:
%   cartesian - A vector or list of vectors [x, y, z] in Cartesian coordinates.
% Output:
%   spherical - A vector or list of vectors [r, theta, phi] in spherical coordinates where:
%               r     - Radial distance (magnitude of the vector).
%               theta - Polar angle (angle from the positive z-axis in radians).
%               phi   - Azimuthal angle (angle from the positive x-axis in radians).

% Convert Cartesian coordinates to spherical coordinates
x = cartesian(:, 1);
y = cartesian(:, 2);
z = cartesian(:, 3);

% Compute radial distance
r = sqrt(x.^2 + y.^2 + z.^2);

% Compute polar angle (theta)
theta = acos(z ./ r);
theta(r == 0) = 0; % Handle the singularity at the origin

% Compute azimuthal angle (phi)
phi = mod(atan2(y, x) + 2*pi, 2*pi);

% Combine into spherical coordinates
spherical = [r, theta, phi];
end