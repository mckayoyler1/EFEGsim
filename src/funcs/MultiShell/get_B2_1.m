% Calculate B2_1 for normalization (equation 16)
function B2_1 = get_B2_1(l, shells)
    % This implements equation 16
    lambda1 = calculate_lambda(l, shells(1).epsilon, shells(1).eta);
    B2_1 = 1 / (shells(1).epsilon * (2*lambda1 - 1));
end