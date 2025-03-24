% Calculate lambda based on equation 10
function lambda = calculate_lambda(l, epsilon, eta)
if epsilon == eta
    lambda = l;
else
    lambda = sqrt(l*(l+1)*eta/epsilon + 1/4) - 1/2;
end
end
