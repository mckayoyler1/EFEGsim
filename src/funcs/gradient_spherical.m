function gradients = gradient_spherical(cartPos, values, distanceCutoff)
% GRADIENT_SPHERICAL_COORDS Computes the gradient of scalar values in spherical coordinates.
%
% Inputs:
%   cartPos       - Nx3 matrix of **Cartesian** coordinates [x, y, z].
%   values        - Nx1 vector of scalar values (e.g., potential values) at the positions.
%   distanceCutoff (optional) - scalar specifying the maximum neighbor distance 
%                               to include in the gradient calculation.
%                               If not provided or set to Inf, no cutoff is applied.
%
% Outputs:
%   gradients - Nx3 matrix of approximate spherical gradients [dr, dtheta, dphi].
%
% NOTE: This is a naive O(N^2) IDW approach in spherical coordinates. For large N,
%       or if you need high accuracy, consider local methods (k-nearest neighbors,
%       local polynomial fits, etc.).

    if nargin < 3 || isempty(distanceCutoff)
        distanceCutoff = Inf;  % Default: no cutoff
    end

    % Validate inputs
    if size(cartPos, 2) ~= 3
        error('Positions must be an Nx3 matrix of Cartesian coordinates [x, y, z].');
    end
    if size(cartPos, 1) ~= numel(values)
        error('Number of positions must match the number of values.');
    end
    
    % Convert from Cartesian to spherical coordinates
    sphPos = toSphere(cartPos);  % [r, theta, phi]
    
    % Extract spherical coordinates
    r     = sphPos(:, 1);     % Radial distances
    theta = sphPos(:, 2);     % Polar angles
    phi   = sphPos(:, 3);     % Azimuthal angles

    % Number of points
    N = size(sphPos, 1);

    % Initialize gradient output
    gradients = zeros(N, 3);

    % Compute the gradient using finite differences (IDW)
    for i = 1:N
        % Current point's spherical coords
        ri     = r(i);
        thetai = theta(i);
        phii   = phi(i);

        % Initialize derivatives
        dr     = 0;
        dtheta = 0;
        dphi   = 0;
        weight_sum = 0;

        for j = 1:N
            if i == j
                continue;  % Skip self
            end
            
            % Neighbor's coords
            rj     = r(j);
            thetaj = theta(j);
            phij   = phi(j);

            % Differences in spherical coords
            delta_r     = rj - ri;
            delta_theta = thetaj - thetai;
            delta_phi   = phij - phii;

            % Handle wrapping of phi (so delta_phi is in [-pi, pi])
            if delta_phi > pi
                delta_phi = delta_phi - 2*pi;
            elseif delta_phi < -pi
                delta_phi = delta_phi + 2*pi;
            end

            % Approximate distance in spherical space
            % (Note this is not true Euclidean distance, but local spherical increments)
            distance = sqrt( delta_r^2 ...
                           + (ri * delta_theta)^2 ...
                           + (ri * sin(thetai) * delta_phi)^2 );

            % Apply the distance cutoff
            if distance < distanceCutoff
                % Inverse-distance-squared weight
                weight = 1 / (distance^2);

                % Increment derivatives based on spherical directions
                df = values(j) - values(i);

                % d/dr
                dr = dr + weight * df * (delta_r / distance);

                % d/dθ (factor of 1/ri)
                if ri > 0
                    dtheta = dtheta + weight * df * (delta_theta / (distance * ri));
                end

                % d/dφ (factor of 1/(ri sin θ))
                if ri > 0 && sin(thetai) ~= 0
                    dphi = dphi + weight * df * (delta_phi / (distance * ri * sin(thetai)));
                end

                % Accumulate total weight
                weight_sum = weight_sum + weight;
            end
        end
        
        % Normalize derivatives by sum of weights
        if weight_sum > 0
            dr     = dr     / weight_sum;
            dtheta = dtheta / weight_sum;
            dphi   = dphi   / weight_sum;
        end

        % Store the gradient for this point
        gradients(i, :) = [dr, dtheta, dphi];
    end
end
