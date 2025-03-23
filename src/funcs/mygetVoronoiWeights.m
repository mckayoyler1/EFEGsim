
function weights = mygetVoronoiWeights(points)
% getVoronoiWeights Compute approximate Voronoi weights for points on the unit sphere.
%
%   weights = getVoronoiWeights(points) returns an Mx1 vector of weights (in steradians)
%   corresponding to the areas of the Voronoi cells on the unit sphere.
%
%   Input:
%       points: an Mx2 matrix of spherical coordinates [theta, phi] in radians,
%               where theta is the polar (inclination) angle (0 to pi) and phi is the
%               azimuth (0 to 2*pi).
%
%   Output:
%       weights: an Mx1 vector containing the area (in steradians) of each Voronoi cell.
%
%   This function works by:
%       1. Converting the spherical coordinates to Cartesian coordinates.
%       2. Computing the convex hull of the points (which yields a triangulation of the sphere).
%       3. Computing the spherical area of each triangle using L'Huilier's formula.
%       4. Accumulating one-third of each triangle's area to each of its vertices.
%
%   Example:
%       % Create a set of random points on the sphere:
%       M = 100;
%       theta = acos(2*rand(M,1)-1);  % polar angle (0 to pi)
%       phi = 2*pi*rand(M,1);         % azimuth (0 to 2*pi)
%       pts = [theta, phi];
%       weights = getVoronoiWeights(pts);
%       disp(weights);
%

    % Convert spherical [theta, phi] to Cartesian coordinates on the unit sphere.
    theta = points(:,1);
    phi = points(:,2);
    x = sin(theta) .* cos(phi);
    y = sin(theta) .* sin(phi);
    z = cos(theta);
    pts_cart = [x, y, z];

    % Compute a triangulation via the convex hull.
    % For points on the sphere in general position, convhulln returns triangles
    % that tessellate the sphere.
    tri = convhulln(pts_cart);
    
    M = size(points, 1);
    weights = zeros(M,1);
    
    % Loop over each triangle in the convex hull.
    % Each row of 'tri' contains indices into pts_cart for a triangle.
    for i = 1:size(tri,1)
        idx = tri(i, :);
        v1 = pts_cart(idx(1), :)';
        v2 = pts_cart(idx(2), :)';
        v3 = pts_cart(idx(3), :)';
        
        % Compute the lengths of the sides (great-circle distances) using the dot product.
        a = acos( max( min(dot(v2, v3), 1), -1) );
        b = acos( max( min(dot(v1, v3), 1), -1) );
        c = acos( max( min(dot(v1, v2), 1), -1) );
        
        s = (a + b + c)/2;  % semi-perimeter
        
        % Use L'Huilier's formula to compute the spherical excess (area of the triangle).
        tanTerm = tan(s/2) * tan((s - a)/2) * tan((s - b)/2) * tan((s - c)/2);
        % Numerical issues might lead to a small negative value; ensure non-negative.
        tanTerm = max(tanTerm, 0);
        E = 4 * atan(sqrt(tanTerm));  % spherical excess (area in steradians)
        
        % Distribute one-third of the area of this triangle to each vertex.
        for j = 1:3
            weights(idx(j)) = weights(idx(j)) + E/3;
        end
    end
    
    % (Optional) Warn if total area is not close to 4*pi (the surface area of a unit sphere).
    totalArea = sum(weights);
    if abs(totalArea - 2*pi) > 1e-2
        warning('Total computed Voronoi area (%.6f) differs from 4*pi.', totalArea);
    end
end