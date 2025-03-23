function data = getSHBasis(data)
% GENERATECOMPLEXSPHERICALHARMONICBASIS Creates a matrix of complex spherical harmonics 
% evaluated at grid points
%
% Inputs:
%   data - Data structure containing:
%          .sphpos - Nx3 array of [r, theta, phi] coordinates
%          .expansionOrder - Maximum degree of spherical harmonics (L)
%
% Outputs:
%   data - Updated structure with added field:
%          .basis - Nx(L+1)^2 complex matrix where each column is a spherical harmonic
%                  evaluated at the grid points, arranged as:
%                  Y_0,0 | Y_1,-1 | Y_1,0 | Y_1,1 | Y_2,-2 | Y_2,-1 | ... | Y_L,L
%
% Note:
%   This function uses the normalized complex spherical harmonics defined as:
%   Y_l,m = sqrt((2l+1)/(4π) * (l-|m|)!/(l+|m|)!) * P_l,|m|(cos θ) * exp(im*φ)
%   where P_l,m are the associated Legendre polynomials.

% Check for required fields
if ~isfield(data, 'sphpos')
    error('Field sphpos not found in data structure');
end

if ~isfield(data, 'expansionOrder')
    error('Field expansionOrder not found in data structure');
end

% Extract spherical coordinates and expansion order
coords = data.sphpos;
L = data.expansionOrder;

% Number of points and total number of spherical harmonics
num_points = size(coords, 1);
num_harmonics = (L + 1)^2;

% Extract theta and phi values
theta = coords(:, 2);  % Polar angle (0 to π)
phi = coords(:, 3);    % Azimuthal angle (0 to 2π)

% Precompute cos(theta) for Legendre polynomial calculation
cos_theta = cos(theta);

% Initialize the basis matrix (complex)
Y = zeros(num_points, num_harmonics);

% Counter for column index
col_idx = 1;

% Loop through each degree l
for l = 0:L
    % Loop through each order m from -l to l
    for m = -l:l
        % Compute associated Legendre polynomial P_l,|m|(cos θ)
        % Note: MATLAB's legendre function uses a different normalization
        % and returns P_l,m for m=0,1,...,l, so we need to handle m<0 separately
        abs_m = abs(m);
        
        % Get the associated Legendre polynomials (MATLAB's normalization)
        P = legendre(l, cos_theta);
        
        % Extract the appropriate order |m|
        if size(P, 1) > 1
            P_lm = P(abs_m+1, :)';  % +1 because MATLAB uses 1-based indexing
        else
            P_lm = P';  % For l=0, m=0
        end
        
        % Apply normalization to convert from MATLAB's convention
        % to normalized complex spherical harmonics
        
        % Apply Condon-Shortley phase
        if m < 0
            cs_phase = (-1)^m;
        else
            cs_phase = 1;
        end
        
        % Normalization factor
        norm_factor = sqrt((2*l+1)/(4*pi) * factorial(l-abs_m)/factorial(l+abs_m));
        
        % Apply the phase and normalization
        P_lm = cs_phase * norm_factor * P_lm;
        
        % Add the complex exponential part
        Y(:, col_idx) = P_lm .* exp(1i * m * phi);
        
        % Move to next column
        col_idx = col_idx + 1;
    end
end

% Store the basis matrix in the data structure
data.basis = Y;

end
