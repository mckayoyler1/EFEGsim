function sh = getSHEntry(SH_list, l, m)
% getSHEntry Extracts the spherical harmonic corresponding to (l, m).
%
%   sh = getSHEntry(SH_list, l, m) returns the spherical harmonic corresponding
%   to degree l and order m from the flattened list SH_list. The ordering is assumed to be:
%
%       [Y_0^0, Y_1^-1, Y_1^0, Y_1^1, Y_2^-2, Y_2^-1, Y_2^0, Y_2^1, Y_2^2, ...].
%
%   The single-index is computed as:
%
%       idx = l*(l+1) + m + 1.
%
%   Inputs:
%       SH_list : Either a vector of length (L+1)^2 (e.g., spherical harmonic coefficients)
%                 or a matrix of size (K x (L+1)^2) where each column corresponds to the
%                 evaluation of one spherical harmonic at K points.
%       l       : Degree (nonnegative integer)
%       m       : Order (integer with -l <= m <= l)
%
%   Output:
%       sh      : The extracted spherical harmonic. If SH_list is a vector, sh is a scalar.
%                 If SH_list is a matrix, sh is a column vector (the corresponding column).
%
%   Example:
%       % Suppose Y_matrix is a 100x9 matrix corresponding to spherical harmonics for l=0,1,2.
%       sh = getSHEntry(Y_matrix, 1, 0);  % Extracts the column for Y_1^0.
%
%       % For a vector of coefficients:
%       coeffs = randn(9,1);
%       a10 = getSHEntry(coeffs, 1, 0);
%

    % Check that m is valid for degree l.
    if m < -l || m > l
        error('Invalid order m: must satisfy -l <= m <= l.');
    end

    % Compute the index using our ordering.
    idx = l*(l+1) + m + 1;
    
    % Extract the entry.
    if isvector(SH_list)
        if idx > length(SH_list)
            error('Index exceeds length of SH_list vector.');
        end
        sh = SH_list(idx);
    else
        % SH_list is a matrix (each column is one harmonic).
        if idx > size(SH_list, 2)
            error('Index exceeds number of columns in SH_list.');
        end
        sh = SH_list(:, idx);
    end
end