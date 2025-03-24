% Get A1, B1 coefficients for a given shell and harmonic index
function [A1, B1] = get_coefficients_A1B1(j, l, shells)
    % This would implement equations 17-20 to calculate A1, B1 recursively
    % Start with innermost shell where A1_1 = 1, B1_1 = 0
    % Then propagate outward using the C matrices
    
    % Initialize with innermost shell
    A1 = 1;
    B1 = 0;
    
    % For shells 2 to j, propagate coefficients
    for i = 1:j-1
        r = shells(i).outer_radius;  % radius at boundary
        lambda_i = calculate_lambda(l, shells(i).epsilon, shells(i).eta);
        lambda_i1 = calculate_lambda(l, shells(i+1).epsilon, shells(i+1).eta);
        
        % Construct C matrices
        C_i = [r^lambda_i, r^(-lambda_i-1); 
               shells(i).epsilon*lambda_i*r^(lambda_i-1), -shells(i).epsilon*(lambda_i+1)*r^(-lambda_i-2)];
        
        C_i1 = [r^lambda_i1, r^(-lambda_i1-1); 
                shells(i+1).epsilon*lambda_i1*r^(lambda_i1-1), -shells(i+1).epsilon*(lambda_i1+1)*r^(-lambda_i1-2)];
        
        % Calculate new coefficients
        coefs = C_i1 \ (C_i * [A1; B1]);
        A1 = coefs(1);
        B1 = coefs(2);
    end
end