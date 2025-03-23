function sphericalHarmonicLSQv2(NTH, JMAX, data, theta, phi)
% Computes spherical harmonic coefficients using least-squares adjustment
%
% Inputs:
%   NTH  - Number of latitude circles (must be even)
%   JMAX - Maximum degree of spherical harmonics
%   data - 2D array of function values on a regular angular grid
%   theta - Vector of theta (colatitude) values [0, pi/2]
%   phi   - Vector of phi (longitude) values [0, 2*pi]
%
% Outputs:
%   Prints spherical harmonic coefficients and reconstructs the data.

    % Validate inputs
    if mod(NTH, 2) ~= 0
        error('NTH (number of latitude circles) must be even.');
    end
    if size(data, 1) ~= length(theta) || size(data, 2) ~= length(phi)
        error('Data dimensions must match length(theta) x length(phi).');
    end

    % Precompute associated Legendre polynomials
    P = computeLegendrePolynomials(JMAX, theta);

    % Compute the right-hand sides of the normal equations
    [CRHS, G] = computeRightHandSidesAndMatrix(theta, phi, JMAX, data, P);

    % Solve the least-squares problem
    COEF = G \ CRHS; % Solve the system with computed normal matrix

    % Synthesize the scalar field from the coefficients
    reconstructedData = synthesizeScalarField(theta, phi, JMAX, COEF, P);

    % Display results
    disp('Spherical harmonic coefficients:');
    for l = 0:JMAX
        for m = 0:l
            idx = l * (l + 1) / 2 + m + 1;
            fprintf('Degree %d, Order %d: %.6e\n', l, m, COEF(idx));
        end
    end

    % Compute reconstruction error
    error = sqrt(mean((data(:) - reconstructedData(:)).^2));
    fprintf('Reconstruction RMS error: %.6e\n', error);
end

function P = computeLegendrePolynomials(JMAX, theta)
% Computes associated Legendre polynomials for all degrees and orders
    x = cos(theta(:)); % Argument for Legendre polynomials
    P = cell(JMAX + 1, JMAX + 1);
    for l = 0:JMAX
        for m = 0:l
            P{l+1, m+1} = legendre(l, x, 'sch');
        end
    end
end

function [CRHS, G] = computeRightHandSidesAndMatrix(theta, phi, JMAX, data, P)
% Computes the right-hand sides and normal matrix for the least-squares problem

    % Initialize CRHS and G
    numCoeffs = (JMAX + 1) * (JMAX + 2) / 2;
    CRHS = zeros(numCoeffs, 1);
    G = zeros(numCoeffs, numCoeffs);

    % Loop over degrees and orders to compute CRHS and G
    for l1 = 0:JMAX
        for m1 = 0:l1
            idx1 = l1 * (l1 + 1) / 2 + m1 + 1;
            for l2 = 0:JMAX
                for m2 = 0:l2
                    idx2 = l2 * (l2 + 1) / 2 + m2 + 1;
                    sumG = 0;
                    for k = 1:length(theta)
                        for j = 1:length(phi)
                            Ylm1 = P{l1+1, m1+1}(k) * exp(1i * m1 * phi(j));
                            Ylm2 = P{l2+1, m2+1}(k) * exp(1i * m2 * phi(j));
                            sumG = sumG + conj(Ylm1) * Ylm2 * sin(theta(k));
                        end
                    end
                    G(idx1, idx2) = sumG;
                end
            end

            % Compute CRHS for the current degree/order
            sumCRHS = 0;
            for k = 1:length(theta)
                for j = 1:length(phi)
                    Ylm = P{l1+1, m1+1}(k) * exp(1i * m1 * phi(j));
                    sumCRHS = sumCRHS + data(k, j) * conj(Ylm) * sin(theta(k));
                end
            end
            CRHS(idx1) = sumCRHS;
        end
    end
end

function reconstructedData = synthesizeScalarField(theta, phi, JMAX, COEF, P)
% Reconstructs the scalar field from spherical harmonic coefficients

    reconstructedData = zeros(length(theta), length(phi));

    for k = 1:length(theta)
        for j = 1:length(phi)
            for l = 0:JMAX
                for m = 0:l
                    idx = l * (l + 1) / 2 + m + 1;
                    Ylm = P{l+1, m+1}(k) * exp(1i * m * phi(j));
                    reconstructedData(k, j) = reconstructedData(k, j) + COEF(idx) * Ylm;
                end
            end
        end
    end
end
