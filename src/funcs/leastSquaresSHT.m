function data = leastSquaresSHT(data)
    
disp('Beginning Least Squares Transform')

    % Check if necessary fields exist
    if ~isfield(data, 'basis')
        builtin('error','Spherical harmonic basis not found. Run generateComplexSphericalHarmonicBasis first.');
    end
    
    % Get the basis and function values
    Y = data.basis;  % Spherical harmonic basis
    YV = data.scaledBasis;
    F = data.V.true;  % Function values at grid points

    % Check if we have weights for weighted least squares
    if isfield(data, 'weights')
        weights = data.weights;

        % % Perform weighted least-squares transform
        % disp('Using weighted least squares with Voronoi weights');
        % coeffs = (Y' * diag(weights) * Y) \ (Y' * diag(weights) * F);

        % Add Tikhonov regularization to least squares
        lambda = 0.0005;  % Regularization parameter
        L = data.expansionOrder;

        % Create regularization matrix (penalizes higher order harmonics more)
        reg_matrix = zeros(size(data.basis, 2));
        idx = 1;
        for l = 0:L
            for m = -l:l
                reg_matrix(idx, idx) = l^2;  % or l^4 for stronger dampening
                idx = idx + 1;
            end
        end

        % Apply regularized least squares
        coeffs = (YV' * diag(weights) *  YV + lambda * reg_matrix) \ (YV' * diag(weights) * F);
    else
        % Perform standard least squares transform
        disp('Using standard least squares (no weights)');
        coeffs = pinv(YV) * F;
    end
    
    % Store the coefficients in the data structure
    data.coeffs = coeffs;
    
    % Calculate goodness of fit
    F_reconstructed = YV * coeffs;
    error = F - F_reconstructed;
    rmse = sqrt(mean(abs(error).^2));
    
    % Store error metrics
    data.reconstructionError = struct('rmse', rmse);
    
    fprintf('Least Squares completed with RMSE: %.6e\n', rmse);
    
    % Return the updated data structure
    return;
end
