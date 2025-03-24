function data = computeSHExpansion(data)
% computeSphericalHarmonicExpansion Compute spherical harmonic expansion info.
%
%   result = computeSphericalHarmonicExpansion(points, V, expansionOrder, basisType)
%
%   INPUTS:
%       points         - N-by-3 array of Cartesian coordinates.
%       V              - N-by-1 vector of measured scalar values.
%       expansionOrder - Order of the spherical harmonic expansion.
%       basisType      - Type of spherical harmonic basis ('complex' or 'real').
%
%   OUTPUT:
%       result         - A struct containing:
%                          .points         - Input Cartesian coordinates.
%                          .sphericalPoints- Spherical coordinates corresponding to points.
%                          .SPHang         - [azimuth, inclination] angles.
%                          .voronoiWeights - Voronoi weights for fitting.
%                          .shExpansion    - Spherical harmonic expansion coefficients.
%                          .reconstructedV - Values reconstructed from the SH expansion.
%                          .expansionOrder - The expansion order used.
%                          .basisType      - The type of spherical harmonic basis used.
%
%   NOTE:
%       This function assumes the existence of helper functions:
%           - toSphere: converts Cartesian to spherical coordinates.
%           - getVoronoiWeights: computes the weights for the fit.
%           - leastSquaresSHT: computes the SH coefficients.
%           - inverseSHT: reconstructs values from the SH coefficients.
%
%   Example:
%       sensorArray = SensorArray(config);
%       points = sensorArray.Locations;
%       V = sensorArray.MeasuredV;
%       result = computeSphericalHarmonicExpansion(points, V, 10, 'complex');
%       % Now result.shExpansion, result.reconstructedV, etc. are available.



    %% Convert Cartesian points to spherical coordinates.
    % Assumption: toSphere returns an N-by-3 array [radius, elevation, azimuth]

    % %% Compute the Voronoi weights.
    % data = calculateVoronoiWeights(data);
    %% Calculate weights
    if strcmp(data.grid_type, 'uniform') || strcmp(data.grid_type, 'modulated')
        data = calculateAnalyticalWeights(data);
    % else
    %     data = calculateVoronoiWeights(data);
    end

    %% Get fields for leastSquares
    Vtrue  = evaluateSymbV(data.pos, data.V.symbV);
    %Vtrue = real(computeYlm(2,2, data.sphpos(:,2), data.sphpos(:,3))/(data.cfg.headshape.radius^(2+1)));
    Etrue  = evaluateSymbE(data.pos, data.E.symbE);
    data.V.true  = Vtrue;
    data.E.true  = Etrue;
    
    data = getScaledSHBasis(data);
    data = leastSquaresSHT(data);
    data = getVfromSHT(data);
    data = getEfromSHT(data);
    
    % %% Package all results into a struct
    % 
    % % data.points         = pos;
    % % data.sphPoints      = sphpos;
    % % data.SPHang         = dirs;
    % % % % data.voronoiWeights = voronoiWeights;
    % % data.coeffs         = coeffs;
    % data.V.true         = Vtrue;
    % data.E.true         = Etrue;
    % % data.expansionOrder = expansionOrder;
    % % data.basisType      = basisType;

end