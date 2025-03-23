function V = computeSymbolicPotential(dipoles, x, y, z)
    % computeSymbolicPotential Symbolically computes the potential at a point (x, y, z)
    % due to a list of dipole objects.
    %
    % Inputs:
    %   dipoles - Array of Dipole objects
    %   x, y, z - Symbolic variables representing the observation point
    %   k       - Constant for potential calculation
    %
    % Output:
    %   V - Symbolic expression for the potential at the observation point

    % Initialize symbolic potential
    V = 0;
    % Loop over all dipoles
    for i = 1:numel(dipoles)
        dipole = dipoles(i);

        % Extract dipole properties
        r_d = dipole.Position; % Dipole position as a vector [x_d, y_d, z_d]
        q_d = dipole.Moment;   % Dipole moment as a vector [qx, qy, qz]

        % Compute symbolic distance vector and magnitude
        r_vec = [x - r_d(1), y - r_d(2), z - r_d(3)]; % Distance vector
        r_mag = sqrt(r_vec(1)^2 + r_vec(2)^2 + r_vec(3)^2); % Distance magnitude

        % Compute symbolic potential due to this dipole
        V = V + Constants.k * dot(q_d, r_vec) / r_mag^3;
    end
end