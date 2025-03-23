function data = getVfromSHT(data)
coeffs = data.coeffs;
a_N = coeffs;
Y_N = data.basis;
Vreconstructed       = Y_N * a_N;
data.V.reconstructed = real(Vreconstructed);
data.V.complexSHT    = Vreconstructed;

end

