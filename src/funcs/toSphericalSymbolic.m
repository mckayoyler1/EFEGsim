function V_spherical = toSphericalSymbolic(V_cartesian)
    % cartesianToSpherical Converts a symbolic function from Cartesian to spherical coordinates
    %
    % Inputs:
    %   V_cartesian - Symbolic function in Cartesian coordinates (x, y, z)
    % Outputs:
    %   V_spherical - Symbolic function in spherical coordinates (r, theta, phi)

    % Define symbolic variables
    syms x y z r theta phi

    % Cartesian to spherical coordinate transformation
    x_sph = r * sin(theta) * cos(phi);
    y_sph = r * sin(theta) * sin(phi);
    z_sph = r * cos(theta);

    % Substitute Cartesian variables with spherical equivalents
    V_spherical = subs(V_cartesian, [x, y, z], [x_sph, y_sph, z_sph]);
end