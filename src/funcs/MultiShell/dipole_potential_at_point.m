function Phi = dipole_potential_at_point(sphpos, r0, alpha, q, model_params)
    % Implement equation 31 from Petrov's paper
    % sphpos = r, theta, phi: observation point coordinates
    % r0: source dipole location (distance from center)
    % alpha: angle between dipole moment q and z axis
    % q: dipole magnitude
    % model_params: structure containing shell parameters
    
    r= sphpos(1);
    theta = sphpos(2);
    phi = sphpos(3);

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
    
    % Calculate potential
    Phi = 0;
    
    % Sum over l (spherical harmonic index)
    for l = 1:50  % Truncate at some reasonable value like 50
        % Calculate lambda for source and observation shells
        lambda_j = calculate_lambda(l, shells(j).epsilon, shells(j).eta);
        lambda_j0 = calculate_lambda(l, shells(j0).epsilon, shells(j0).eta);
        
        % Get A1, B1 for source shell and A2, B2 for observation shell
        [A1_j0, B1_j0] = get_coefficients_A1B1(j0, l, shells);
        [A2_j, B2_j] = get_coefficients_A2B2(j, l, shells);
        
        % Calculate B2_1 (for innermost shell)
        B2_1 = get_B2_1(l, shells);
        
        % Calculate factor in front of the sum
        factor = 1 / (4 * pi * epsilon1) * (1 / B2_1) * ((2*l + 1) / (2*lambda_j0 + 1));
        
        % Calculate potential contribution
        term1 = (A2_j * r^lambda_j + B2_j / r^(lambda_j+1));
        
        % qr term (radial dipole component)
        qr_term = qr * (lambda_j0 * A1_j0 * r0^(lambda_j0-1) - (lambda_j0 + 1) * B1_j0 / r0^(lambda_j0+2)) * legendre_P(l, 0, cos(theta));
        
        % qtheta term (tangential dipole component)
        qtheta_term = -qtheta * (A1_j0 * r0^(lambda_j0-1) + B1_j0 / r0^(lambda_j0+2)) * legendre_P(l, 1, cos(theta)) * cos(phi);
        
        % Add this harmonic's contribution
        Phi = Phi + factor * term1 * (qr_term + qtheta_term);
    end
end

