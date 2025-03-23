function points = unflattenMeshgrid(flatPoints, vals, numTheta, numPhi)
    % unflattenMeshgrid converts a flattened Nx3 array back into meshgrid form.
    %
    % Inputs:
    %   flatPoints - Nx3 array of points [x, y, z].
    %
    % Outputs:
    %   X, Y, Z    - 3D arrays corresponding to the original meshgrid.

    % Validate inputs
    if size(flatPoints, 2) ~= 3
        error('flatPoints must be an Nx3 array.');
    end

    flatPoints = toSphere(flatPoints);
    flatPoints = [flatPoints vals];
    points = reshape(flatPoints,numTheta,numPhi,4);
    
    
end