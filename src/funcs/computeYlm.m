function Ylm = computeYlm(l, m, theta, phi)
% computeYlmNumeric Computes the complex spherical harmonic Y_l^m at given angles.
%
%   Ylm = computeYlmNumeric(l, m, theta, phi) returns a column vector of values
%   of the complex spherical harmonic of degree l and order m evaluated at the angles
%   theta and phi. theta and phi must be vectors of the same length.
%
%   The spherical harmonic is defined as:
%
%     Y_l^m(theta,phi) = sqrt((2*l+1)/(4*pi) * (l-|m|)!/(l+|m|)!) * P_l^{|m|}(cos(theta)) * exp(1i*m*phi)
%
%   and for negative m, the identity
%
%     Y_l^{-m} = (-1)^m * conj(Y_l^m)
%
%   is used.

    mabs = abs(m);
    % Compute the associated Legendre polynomials.
    % MATLAB's legendre returns an (l+1) x numTheta array.
    x = cos(theta);  % vector of cos(theta)
    Pl_all = legendre(l, x);  % size: (l+1) x length(x)
    % Get the row corresponding to order mabs.
    Plm = Pl_all(mabs+1, :)';
    
    % Normalization factor:
    normFactor = sqrt((2*l+1)/(4*pi) * factorial(l-mabs)/factorial(l+mabs));
    
    if m >= 0
        Ylm = normFactor .* Plm .* exp(1i*m*phi);
    else
        % For negative m, use Y_l^{-m} = (-1)^m * conj(Y_l^{|m|})
        Ylm = ((-1)^m) .* conj(normFactor .* Plm .* exp(1i*mabs*phi));
    end
end
