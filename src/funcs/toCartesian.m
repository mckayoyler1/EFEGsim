function [x,y,z] = toCartesian(spherical)
            % Convert spherical coordinates to Cartesian coordinates
            
            r = spherical(:, 1);
            theta = spherical(:, 2);
            phi = spherical(:, 3);
            rsin = r.* sin(theta);
            x = rsin .* cos(phi);
            y = rsin .* sin(phi);
            z = r .* cos(theta);
end