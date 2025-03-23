function F = inverseSHT(data)
%INVERSE_SHT Perform the inverse spherical harmonic transform
%
%   F_N:    (N+1)^2 x L matrix of SH coefficients up to order N,  
%           with L spherical functions encoded as columns


    Y_N = data.basis;
    F_N = data.coeffs;
    % perform the inverse transform up to degree N
    F = Y_N * F_N;

end
