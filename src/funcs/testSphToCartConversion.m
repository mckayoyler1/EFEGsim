function testSphToCartConversion(tol)
% checkSphToCartConversion Checks the spherical-to-cartesian vector field conversion.
%
%   checkSphToCartConversion(tol) creates a sample spherical vector field (purely
%   radial) and compares the converted Cartesian vectors to the expected values.
%
%   tol is an optional tolerance (default is 1e-10).
%
%   The function uses convertSphericalToCartesian to perform the conversion.
%
%   Example:
%       checkSphToCartConversion(1e-10);

    if nargin < 1
        tol = 1e-10;
    end

    % Create a grid of spherical points on the unit sphere.
    [phi, theta] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 50));
    r = ones(size(theta));  % Unit sphere

    % Reshape the spherical points into an Mx3 matrix: [r, theta, phi]
    sphPoints = [r(:), theta(:), phi(:)];

    % Define a purely radial vector field: Vr = 1, Vtheta = 0, Vphi = 0.
    V_r     = ones(size(r(:)));
    V_theta = zeros(size(r(:)));
    V_phi   = zeros(size(r(:)));
    sphVectors = [V_r, V_theta, V_phi];

    % Convert the spherical vector field to Cartesian coordinates.
    cartVectors = convertSphericalVecToCartesian(sphPoints, sphVectors);

    % Expected Cartesian vector for a purely radial field:
    % [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    expectedVx = sin(sphPoints(:,2)) .* cos(sphPoints(:,3));
    expectedVy = sin(sphPoints(:,2)) .* sin(sphPoints(:,3));
    expectedVz = cos(sphPoints(:,2));
    expectedCartVectors = [expectedVx, expectedVy, expectedVz];

    % Calculate the maximum absolute error.
    maxError = max(abs(cartVectors(:) - expectedCartVectors(:)));

    % Check the error against the tolerance.
    if maxError < tol
        fprintf('Conversion check passed. Maximum error: %e\n', maxError);
    else
        error('Conversion check failed. Maximum error: %e\n', maxError);
    end
end
