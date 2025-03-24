function Efield = dipole_electric_field_at_point(sphpoint, r0, alpha, q, model_params)
    % Calculates electric field components at a given point
    % r, theta, phi: observation point coordinates
    % r0: source dipole location (distance from center)
    % alpha: angle between dipole moment q and z axis
    % q: dipole magnitude
    % model_params: structure containing shell parameters
    r = sphpoint(1);
    theta = sphpoint(2);
    phi = sphpoint(3);

    % Extract model parameters
    epsilon1 = model_params.epsilon1;  % conductivity of innermost shell
    shells = model_params.shells;      % shell parameters
    
    % Calculate dipole components
    qr = q * cos(alpha);  % radial component
    qtheta = q * sin(alpha);  % tangential component
    
    % Find which shell contains the observation point
    j = find_shell_index(r, shells);
    
    % Find which shell contains the source
    j0 = find_shell_index(r0, shells);
    
    % Initialize field components
    Er = 0;
    Etheta = 0;
    Ephi = 0;
    
    % Sum over l (spherical harmonic index) - start from l=2
    for l = 2:50  % Truncate at some reasonable value like 50
        % Calculate lambda for source and observation shells
        lambda_j = calculate_lambda(l, shells(j).epsilon, shells(j).eta);
        lambda_j0 = calculate_lambda(l, shells(j0).epsilon, shells(j0).eta);
        
        % Get coefficients for source shell and observation shell
        [A1_j0, B1_j0] = get_coefficients_A1B1(j0, l, shells);
        [A2_j, B2_j] = get_coefficients_A2B2(j, l, shells);
        
        % Calculate B2_1 (for innermost shell)
        B2_1 = get_B2_1(l, shells);
        
        % Common factor for all field components
        common_factor = 1 / (4 * pi * epsilon1) * (1 / B2_1) * ((2*l + 1) / (2*lambda_j0 + 1));
        
        % Source terms
        source_qr_term = qr * (lambda_j0 * A1_j0 * r0^(lambda_j0-1) - (lambda_j0 + 1) * B1_j0 / r0^(lambda_j0+2));
        source_qtheta_term = qtheta * (A1_j0 * r0^(lambda_j0-1) + B1_j0 / r0^(lambda_j0+2));
        
        % === Radial component (Er) - Equation 32 ===
        % Field factor for radial component
        Er_field_factor = lambda_j * A2_j * r^(lambda_j-1) - (lambda_j + 1) * B2_j / r^(lambda_j+2);
        
        % Add contribution to Er
        Er = Er - common_factor * Er_field_factor * (...
           source_qr_term * legendre_P(l, 0, cos(theta)) - ... 
            source_qtheta_term * legendre_P(l, 1, cos(theta)) * cos(phi)...
        );
        
        % === Theta component (Etheta) - Equation 33 ===
        % Field factor for theta component
        Etheta_field_factor = A2_j * r^(lambda_j-1) + B2_j / r^(lambda_j+2);
        
        % Add contribution to Etheta
        Etheta = Etheta - common_factor * Etheta_field_factor * (...
            source_qr_term * legendre_P(l, 1, cos(theta)) - ...
            source_qtheta_term * 0.5 * (legendre_P(l, 2, cos(theta)) - l*(l+1)*legendre_P(l, 0, cos(theta))) * cos(phi)...
        );
        
        % === Phi component (Ephi) - Equation 34 ===
        % Only depends on qtheta
        Ephi = Ephi - common_factor * Etheta_field_factor * ... 
            source_qtheta_term * legendre_P(l, 1, cos(theta)) / sin(theta) * sin(phi);
    end
    Efield = [Er, Etheta, Ephi];
end