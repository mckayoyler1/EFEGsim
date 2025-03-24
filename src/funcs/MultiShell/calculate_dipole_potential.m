function Phi_values = calculate_dipole_potential(sphpos, r0, alpha, q, model_params)
    
    r_list = sphpos(:,1);
    theta_list = sphpos(:,2);
    phi_list = sphpos(:,3);

    % Get number of points
    num_points = length(r_list);
    Phi_values = zeros(size(r_list));
    
    % Extract model parameters
    epsilon1 = model_params.epsilon1;
    shells = model_params.shells;
    
    % Calculate dipole components
    qr = q * cos(alpha);
    qtheta = q * sin(alpha);
    
    % Find which shell contains the source
    j0 = find_shell_index(r0, shells);
    
    % Maximum l for spherical harmonic summation
    l_max = 50;
    
    % Precompute source-dependent terms
    lambda_j0_values = zeros(l_max, 1);
    A1_j0_values = zeros(l_max, 1);
    B1_j0_values = zeros(l_max, 1);
    B2_1_values = zeros(l_max, 1);
    
    for l = 1:l_max
        lambda_j0_values(l) = calculate_lambda(l, shells(j0).epsilon, shells(j0).eta);
        [A1_j0_values(l), B1_j0_values(l)] = get_coefficients_A1B1(j0, l, shells);
        B2_1_values(l) = get_B2_1(l, shells);
    end
    
    % Process each point
    for i = 1:num_points
        r = r_list(i);
        theta = theta_list(i);
        phi = phi_list(i);
        
        % Find which shell contains the observation point
        j = find_shell_index(r, shells);
        
        % Calculate potential for this point
        Phi = 0;
        
        for l = 1:l_max
            lambda_j = calculate_lambda(l, shells(j).epsilon, shells(j).eta);
            [A2_j, B2_j] = get_coefficients_A2B2(j, l, shells);
            
            % Calculate factor for this harmonic
            factor = 1 / (4 * pi * epsilon1) * (1 / B2_1_values(l)) * ((2*l + 1) / (2*lambda_j0_values(l) + 1));
            
            % Calculate observation point terms
            term1 = (A2_j * r^lambda_j + B2_j / r^(lambda_j+1));
            
            % Calculate source terms
            lambda_j0 = lambda_j0_values(l);
            A1_j0 = A1_j0_values(l);
            B1_j0 = B1_j0_values(l);
            
            % qr term (radial dipole component)
            qr_term = qr * (lambda_j0 * A1_j0 * r0^(lambda_j0-1) - (lambda_j0 + 1) * B1_j0 / r0^(lambda_j0+2)) * legendre_P(l, 0, cos(theta));
            
            % qtheta term (tangential dipole component)
            qtheta_term = -qtheta * (A1_j0 * r0^(lambda_j0-1) + B1_j0 / r0^(lambda_j0+2)) * legendre_P(l, 1, cos(theta)) * cos(phi);
            
            % Add this harmonic's contribution
            Phi = Phi + factor * term1 * (qr_term + qtheta_term);
        end
        
        Phi_values(i) = Phi;
    end
end

