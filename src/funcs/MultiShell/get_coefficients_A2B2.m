% Get A2, B2 coefficients for a given shell and harmonic index
function [A2, B2] = get_coefficients_A2B2(j, l, shells)
    % Similar to A1B1 but starts with outermost shell where A2_N = 0, B2_N = 1
    % Then propagate inward
    
    N = length(shells);  % Number of shells
    
    % Initialize with outermost shell
    A2 = 0;
    B2 = 1;
    
    % For shells N-1 down to j, propagate coefficients
    for i = N:-1:j+1
        r = shells(i-1).outer_radius;  % radius at boundary
        lambda_im1 = calculate_lambda(l, shells(i-1).epsilon, shells(i-1).eta);
        lambda_i = calculate_lambda(l, shells(i).epsilon, shells(i).eta);
        
        % Construct C matrices
        C_im1 = [r^lambda_im1, r^(-lambda_im1-1); 
                 shells(i-1).epsilon*lambda_im1*r^(lambda_im1-1), -shells(i-1).epsilon*(lambda_im1+1)*r^(-lambda_im1-2)];
        
        C_i = [r^lambda_i, r^(-lambda_i-1); 
               shells(i).epsilon*lambda_i*r^(lambda_i-1), -shells(i).epsilon*(lambda_i+1)*r^(-lambda_i-2)];
        
        % Calculate new coefficients
        coefs = C_im1 \ (C_i * [A2; B2]);
        A2 = coefs(1);
        B2 = coefs(2);
    end
end