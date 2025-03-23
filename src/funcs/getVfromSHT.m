function data = getVfromSHT(data)
coeffs = data.coeffs;
r = data.sphpos(:,1);
L = data.expansionOrder;
V_rescaled_coeffs = zeros(size(coeffs));
E_rescaled_coeffs = zeros(size(coeffs));
for l = 0:L
    for m = -l:l
        idx = l^2 + l + m + 1;  % Index in the coefficient vector
        V_rescaled_coeffs(idx) = coeffs(idx) ./ (r(idx)^(l+1));
        E_rescaled_coeffs(idx) = coeffs(idx) ./ (r(idx)^(l+2));
    end
end
a_N = V_rescaled_coeffs;
Y_N = data.basis;

Vreconstructed       = Y_N * a_N;
data.V.reconstructed = real(Vreconstructed);
data.V.complexSHT    = Vreconstructed;
data.V_ScaledCoeffs  = V_rescaled_coeffs;
data.E_ScaledCoeffs  = E_rescaled_coeffs;
end

