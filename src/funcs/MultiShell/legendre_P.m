% Legendre polynomial function 
function P = legendre_P(l, m, x)
    % For m=0 and m=1, extract the appropriate associated Legendre polynomial
    P_all = legendre(l, x);
    P = P_all(m+1, :);  % m+1 because MATLAB uses 1-based indexing
end