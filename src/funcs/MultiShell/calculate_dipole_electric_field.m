function Efield = calculate_dipole_electric_field(sphpos, r0, alpha, q, model_params)
    r_list = sphpos(:,1);
    theta_list = sphpos(:,2);
    phi_list = sphpos(:,3);
    % Check that all coordinate arrays have the same length
    if length(r_list) ~= length(theta_list) || length(r_list) ~= length(phi_list)
        error('Coordinate arrays must have the same length');
    end
    
    % Initialize results arrays
    num_points = length(r_list);
    Er = zeros(size(r_list));
    Etheta = zeros(size(r_list));
    Ephi = zeros(size(r_list));
    
    % Calculate field components at each point
    for i = 1:num_points
        Efield_i = dipole_electric_field_at_point(sphpos(i,:), r0, alpha, q, model_params);
        Er(i) = Efield_i(1);
        Etheta(i) = Efield_i(2);
        Ephi(i) = Efield_i(3);
    end
    Ephi = replaceNaNwithZero(Ephi);
    Efield = [Er(:), Etheta(:), Ephi(:)];
end
