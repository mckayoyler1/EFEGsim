function [gx,gy,gz] = toCartesianGradient(gradients, cartPos)
    % Convert positions to spherical coordinates
    sphPositions = toSphere(cartPos);
    r = sphPositions(:, 1);
    theta = sphPositions(:, 2);
    phi = sphPositions(:, 3);

    % Convert gradients from spherical to Cartesian
    % dr is radial, dtheta is polar, dphi is azimuthal
    gx = gradients(:, 1) .* sin(theta) .* cos(phi) + ...
         (1 ./ r) .* gradients(:, 2) .* cos(theta) .* cos(phi) - ...
         (1 ./ (r .* sin(theta))) .* gradients(:, 3) .* sin(phi);

    gy = gradients(:, 1) .* sin(theta) .* sin(phi) + ...
         (1 ./ r) .* gradients(:, 2) .* cos(theta) .* sin(phi) + ...
         (1 ./ (r .* sin(theta))) .* gradients(:, 3) .* cos(phi);

    gz = gradients(:, 1) .* cos(theta) - ...
         (1 ./ r) .* gradients(:, 2) .* sin(theta);

end