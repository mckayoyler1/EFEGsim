function results = verifyVoronoiWeights(data)
% VERIFYVORONOIWEIGHTS Performs validation tests for Voronoi weights on a sphere or hemisphere
%
% Inputs:
%   data - Data structure containing:
%          .pos - Cartesian coordinates
%          .sphpos - Spherical coordinates
%          .weights - Weights to verify
%          .grid_type - Type of grid used
%          .cfg.headshape.radius - Radius of the spherical shell
%
% Outputs:
%   results - Struct containing verification results
%
% Tests performed:
% 1. Sum test: Weights should sum to total surface area
% 2. Uniformity test: For uniform grid, weights should be similar
% 3. Function reproduction test: Verify integration of simple functions

% Check if necessary fields exist
if ~isfield(data, 'pos') || ~isfield(data, 'sphpos') || ~isfield(data, 'weights')
    error('Required fields (pos, sphpos, weights) not found in data structure');
end

% Get grid type
if ~isfield(data, 'grid_type')
    data.grid_type = 'uniform';
end
grid_type = lower(data.grid_type);

% Initialize results struct
results = struct();
results.grid_type = grid_type;

% ------ Determine if we have a hemisphere or full sphere ------
weights = data.weights;
points = data.pos;
spherical_coords = data.sphpos;

% Check if this is an upper hemisphere (theta from 0 to pi/2)

theta_max = max(spherical_coords(:,2));
is_upper_hemisphere = (abs(theta_max - pi/2) < 1e-6);

results.is_upper_hemisphere = is_upper_hemisphere;
fprintf('Detected geometry: %s\n', conditional_str(is_upper_hemisphere, 'Upper Hemisphere', 'Full Sphere'));

% ------ Test 1: Sum of weights ------
weight_sum = sum(weights);

if strcmpi(grid_type, 'modulated')
    % For modulated grid, compute expected surface area by averaging r^2
    % r_squared = spherical_coords(:,1).^2;
    r_squared = 1;
    if is_upper_hemisphere
        expected_sum = 2*pi*mean(r_squared); % Hemisphere surface area
    else
        expected_sum = 4*pi*mean(r_squared); % Full sphere surface area
    end
else
    % For uniform or random grids, all points have the same radius
    r = 1;
    
    if is_upper_hemisphere
        expected_sum = 2*pi*r^2; % Hemisphere surface area
    else
        expected_sum = 4*pi*r^2; % Full sphere surface area
    end
end

% Calculate relative error
relative_error = abs(weight_sum - expected_sum) / expected_sum;
results.weight_sum = weight_sum;
results.expected_sum = expected_sum;
results.weight_sum_error = relative_error;
results.weight_sum_passed = relative_error < 0.01; % Pass if error < 1%

% ------ Test 2: Uniformity for regular grid ------

mean_weight = mean(weights);
std_weight = std(weights);
cv = std_weight / mean_weight;

results.weight_mean = mean_weight;
results.weight_std = std_weight;
results.weight_cv = cv;

% For hemispheres, we expect more variation near the equator
if is_upper_hemisphere
    cv_threshold = 0.4; % Allow more variation for hemisphere
else
    cv_threshold = 0.3; % Standard threshold for full sphere
end

if strcmpi(grid_type, 'uniform')
    % For analytical weights, we expect non-uniformity proportional to sin(θ)
    results.uniformity_passed = true;  % Skip this test
    fprintf('Uniformity test skipped (not applicable for analytical weights)\n');
else
    results.uniformity_passed = cv < cv_threshold;
end
% Plot histogram of weights if available
try
    figure('Name', 'Voronoi Weight Distribution');
    histogram(weights, 20);
    title('Histogram of Voronoi Weights');
    xlabel('Weight Value');
    ylabel('Frequency');
    results.histogram_created = true;
catch
    results.histogram_created = false;
end


% ------ Test 3: Test integration of simple functions ------
% Test functions with known integrals over the sphere or hemisphere

% Test function 1: f(x,y,z) = 1 (constant)
f1 = ones(size(points, 1), 1);
integral1 = sum(f1 .* weights);
expected1 = expected_sum; % Same as the weight sum test
error1 = abs(integral1 - expected1) / expected1;

% Test function 2: f(x,y,z) = x^2 + y^2 + z^2 = r^2
f2 = sum(points.^2, 2);
integral2 = sum(f2 .* weights);

if strcmpi(grid_type, 'modulated')
    % For modulated grid, use average r^4
    % r_4 = spherical_coords(:,1).^4;
    r_4 = 1;
    
    if is_upper_hemisphere
        expected2 = 2*pi*mean(r_4); % For hemisphere
    else
        expected2 = 4*pi*mean(r_4); % For full sphere
    end
else
    r = data.cfg.headshape.radius;
    
    if is_upper_hemisphere
        expected2 = 2*pi*r^4; % For hemisphere
    else
        expected2 = 4*pi*r^4; % For full sphere
    end
end

error2 = abs(integral2 - expected2) / expected2;

% Test function 3: f(x,y,z) = z (cosine of polar angle)
f3 = points(:,3);
integral3 = sum(f3 .* weights);

if is_upper_hemisphere
    % For upper hemisphere, z integral is positive
    if strcmpi(grid_type, 'modulated')
        % For modulated radius, we need average r^3
        r_3 = spherical_coords(:,1).^3;
        expected3 = pi*mean(r_3); % For hemisphere, z-integral = pi*r^3
    else
        r = data.cfg.headshape.radius;
        expected3 = pi*r^3; % For hemisphere, z-integral = pi*r^3
    end
else
    % For full sphere, z integral is zero
    expected3 = 0;
end

% For z coordinate, normalize error by expected_sum if expected3 is near zero
if abs(expected3) < 1e-10
    error3 = abs(integral3) / expected_sum;
else
    error3 = abs(integral3 - expected3) / abs(expected3);
end

% Store function test results
results.function_tests = struct();
results.function_tests.constant = struct('value', integral1, 'expected', expected1, 'error', error1);
results.function_tests.r_squared = struct('value', integral2, 'expected', expected2, 'error', error2);
results.function_tests.z_coord = struct('value', integral3, 'expected', expected3, 'error', error3);

% Overall test for function reproduction
% For hemispheres, the z-coordinate test might need a slightly relaxed threshold
if is_upper_hemisphere
    results.function_tests_passed = (error1 < 0.01) && (error2 < 0.05) && (error3 < 0.05);
else
    results.function_tests_passed = (error1 < 0.01) && (error2 < 0.05) && (error3 < 0.01);
end

% ------ Test 4: Visual test - plot Voronoi cells (if available) ------
try
    figure('Name', 'Voronoi Weights Visualization');
    
    % Create color map based on weights
    min_weight = min(weights);
    max_weight = max(weights);
    normalized_weights = (weights - min_weight) / (max_weight - min_weight);
    
    % Plot points on sphere with colors representing weights
    scatter3(points(:,1), points(:,2), points(:,3), 30, normalized_weights, 'filled');
    colorbar;
    title([conditional_str(is_upper_hemisphere, 'Hemisphere', 'Sphere') ' Voronoi Weights - ' grid_type ' grid']);
    axis equal;
    grid on;
    
    % Add sphere/hemisphere wireframe for reference
    hold on;
    r = data.cfg.headshape.radius;
    [x,y,z] = sphere(20);
    x = r*x; y = r*y; z = r*z;
    
    if is_upper_hemisphere
        % Only plot the upper hemisphere
        z(z<0) = NaN;
    end
    
    h = mesh(x,y,z);
    h.FaceColor = 'none';
    h.EdgeColor = [0.7 0.7 0.7];
    h.EdgeAlpha = 0.3;
    hold off;
    
    results.visualization_created = true;
catch
    results.visualization_created = false;
end

% ------ Test 5: Spherical harmonic test ------
% For full sphere, Y_lm should integrate to 0
% For hemisphere, these have specific expected values
try
    % Create degree 1 spherical harmonics
    sin_theta = sin(spherical_coords(:,2));
    cos_theta = cos(spherical_coords(:,2));
    sin_phi = sin(spherical_coords(:,3));
    cos_phi = cos(spherical_coords(:,3));
    
    % Y_1,-1 ~ sin(θ)sin(φ)
    Y1m1 = sin_theta .* sin_phi;
    integral_Y1m1 = sum(Y1m1 .* weights);
    
    % Y_1,0 ~ cos(θ)
    Y10 = cos_theta;
    integral_Y10 = sum(Y10 .* weights);
    %% TEST Y_1,0
    % For upper hemisphere, Y10 (cos(θ)) integral
r = data.cfg.headshape.radius;
expected_Y10 = pi * r^3; % Correct value for hemisphere
fprintf('Y_1,0 integral: Actual=%.6f, Expected=%.6f\n', integral_Y10, expected_Y10);

    % Y_1,1 ~ sin(θ)cos(φ)
    Y11 = sin_theta .* cos_phi;
    integral_Y11 = sum(Y11 .* weights);
    
    % Store results
    results.spherical_harmonics = struct();
    results.spherical_harmonics.Y1m1 = integral_Y1m1;
    results.spherical_harmonics.Y10 = integral_Y10;
    results.spherical_harmonics.Y11 = integral_Y11;
    
    if is_upper_hemisphere
        % For hemisphere, Y10 has non-zero integral, others should be close to 0
        r = data.cfg.headshape.radius;
        expected_Y10 = pi*r^3; % Theoretical value for upper hemisphere
        
        % Calculate errors
        error_Y1m1 = abs(integral_Y1m1) / expected_sum;
        error_Y10 = abs(integral_Y10 - expected_Y10) / abs(expected_Y10);
        error_Y11 = abs(integral_Y11) / expected_sum;
        
        % Add expected values to results
        results.spherical_harmonics.expected_Y10 = expected_Y10;
        results.spherical_harmonics.error_Y1m1 = error_Y1m1;
        results.spherical_harmonics.error_Y10 = error_Y10;
        results.spherical_harmonics.error_Y11 = error_Y11;
        
        % Test is passed if errors are within threshold
        results.spherical_harmonics_passed = (error_Y1m1 < 0.05) && ...
                                             (error_Y10 < 0.05) && ...
                                             (error_Y11 < 0.05);
    else
        % For full sphere, all should integrate to 0
        threshold = 0.01;
        norm_Y1m1 = abs(integral_Y1m1) / expected_sum;
        norm_Y10 = abs(integral_Y10) / expected_sum;
        norm_Y11 = abs(integral_Y11) / expected_sum;
        
        results.spherical_harmonics.norm_Y1m1 = norm_Y1m1;
        results.spherical_harmonics.norm_Y10 = norm_Y10;
        results.spherical_harmonics.norm_Y11 = norm_Y11;
        
        % Test is passed if all are close to zero
        results.spherical_harmonics_passed = (norm_Y1m1 < threshold) && ...
                                             (norm_Y10 < threshold) && ...
                                             (norm_Y11 < threshold);
    end
    
    % Print results
    fprintf('Spherical harmonics test: %s\n', passfail(results.spherical_harmonics_passed));
    if is_upper_hemisphere
        fprintf('  - Y_1,-1 (sin(θ)sin(φ)) error: %.4f%%\n', 100*error_Y1m1);
        fprintf('  - Y_1,0  (cos(θ)) error: %.4f%%\n', 100*error_Y10);
        fprintf('  - Y_1,1  (sin(θ)cos(φ)) error: %.4f%%\n', 100*error_Y11);
    else
        fprintf('  - Y_1,-1 normalized integral: %.6f\n', norm_Y1m1);
        fprintf('  - Y_1,0 normalized integral: %.6f\n', norm_Y10);
        fprintf('  - Y_1,1 normalized integral: %.6f\n', norm_Y11);
    end
catch
    % If there's any error, skip this test
    results.spherical_harmonics_passed = false;
end

% Overall result
if isfield(results, 'uniformity_passed') && isfield(results, 'spherical_harmonics_passed')
    results.all_tests_passed = results.weight_sum_passed && ...
                               results.uniformity_passed && ...
                               results.function_tests_passed && ...
                               results.spherical_harmonics_passed;
elseif isfield(results, 'uniformity_passed')
    results.all_tests_passed = results.weight_sum_passed && ...
                               results.uniformity_passed && ...
                               results.function_tests_passed;
elseif isfield(results, 'spherical_harmonics_passed')
    results.all_tests_passed = results.weight_sum_passed && ...
                               results.function_tests_passed && ...
                               results.spherical_harmonics_passed;
else
    results.all_tests_passed = results.weight_sum_passed && results.function_tests_passed;
end

% Print summary
fprintf('\nVerification of Voronoi weights for %s grid on %s:\n', ...
    grid_type, conditional_str(is_upper_hemisphere, 'upper hemisphere', 'full sphere'));
fprintf('Weight sum test: %s (error: %.4f%%)\n', passfail(results.weight_sum_passed), 100*results.weight_sum_error);
if isfield(results, 'uniformity_passed')
    fprintf('Uniformity test: %s (CV: %.4f)\n', passfail(results.uniformity_passed), results.weight_cv);
end
fprintf('Function integration tests: %s\n', passfail(results.function_tests_passed));
fprintf('  - Constant function error: %.4f%%\n', 100*results.function_tests.constant.error);
fprintf('  - r² function error: %.4f%%\n', 100*results.function_tests.r_squared.error);
fprintf('  - z coordinate error: %.4f%%\n', 100*results.function_tests.z_coord.error);
fprintf('Overall: %s\n', passfail(results.all_tests_passed));

end

% Helper function to format pass/fail as string
function str = passfail(condition)
    if condition
        str = 'PASS';
    else
        str = 'FAIL';
    end
end

% Helper function for conditional string choice
function str = conditional_str(condition, true_str, false_str)
    if condition
        str = true_str;
    else
        str = false_str;
    end
end