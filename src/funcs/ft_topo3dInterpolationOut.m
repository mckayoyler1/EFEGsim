function tri = ft_topo3dInterpolationOut(pos, val, varargin)

% FT_PLOT_TOPO3D visualizes a 3D topographic representation of the electric potential
% or magnetic field distribution at the sensor locations.
%
% Use as
%   ft_plot_topo3d(pos, val, ...)
% where the channel positions are given as a Nx3 matrix and the values are
% given as Nx1 vector.
%
% Optional input arguments should be specified in key-value pairs and can include
%   'contourstyle'  = string, 'none', 'black', 'color' (default = 'none')
%   'isolines'      = vector with values at which to draw isocontours, or 'auto' (default = 'auto')
%   'facealpha'     = scalar, between 0 and 1 (default = 1)
%   'refine'        = scalar, number of refinement steps for the triangulation, to get a smoother interpolation (default = 0)
%   'neighbourdist' = number, maximum distance between neighbouring sensors (default is automatic)
%   'unit'          = string, 'm', 'cm' or 'mm' (default = 'cm')
%   'coordsys'      = string, assume the data to be in the specified coordinate system (default = 'unknown')
%   'axes'          = boolean, whether to plot the axes of the 3D coordinate system (default = false)
%
% See also FT_PLOT_TOPO, FT_PLOT_SENS, FT_PLOT_MESH, FT_PLOT_HEADSHAPE,
% FT_TOPOPLOTER, FT_TOPOPLOTTFR

% get the optional input arguments
contourstyle  = ft_getopt(varargin, 'contourstyle', 'none');
refine_       = ft_getopt(varargin, 'refine', 0);           % do not confuse with the REFINE function
neighbourdist = ft_getopt(varargin, 'neighbourdist');
isolines      = ft_getopt(varargin, 'isolines', 'auto');
topostyle     = ft_getopt(varargin, 'topostyle', 'color');  % FIXME what is the purpose of this option?
facealpha     = ft_getopt(varargin, 'facealpha', 1);
unit          = ft_getopt(varargin, 'unit', 'cm');
coordsys      = ft_getopt(varargin, 'coordsys');
axes_         = ft_getopt(varargin, 'axes', false);         % do not confuse with built-in function

if islogical(contourstyle) && contourstyle==false
  % false was supported up to 18 November 2013, 'none' is more consistent with other plotting options
  contourstyle = 'none';
end

% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

if size(val,2)==size(pos,1)
  val = val';
end

% the interpolation requires a triangulation
tri = projecttri(pos, 'delaunay');

if isempty(neighbourdist)
  % compute the distance between sensor locations
  neighbourdist = min(dist(pos')+diag(inf(size(pos,1),1)), [], 2);
  neighbourdist = 2*max(neighbourdist);
end

if neighbourdist>0 && neighbourdist<inf
  % compute the length of the triangle edges
  v1 = tri(:,1);
  v2 = tri(:,2);
  v3 = tri(:,3);
  len1 = sqrt(sum((pos(v1,:)-pos(v2,:)).^2, 2));
  len2 = sqrt(sum((pos(v2,:)-pos(v3,:)).^2, 2));
  len3 = sqrt(sum((pos(v3,:)-pos(v1,:)).^2, 2));

  % remove triangles with edges that are too long
  skip = any([len1 len2 len3]>neighbourdist, 2);
  tri = tri(~skip,:);
end


if ~isequal(topostyle, false)
  switch topostyle
    case 'color'
      % plot a 2D or 3D triangulated surface with linear interpolation
      if length(val)==size(pos,1)
        hs = patch('Vertices', pos, 'Faces', tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
      else
        hs = patch('Vertices', pos, 'Faces', tri, 'CData', val, 'FaceColor', 'flat');
      end
      set(hs, 'EdgeColor', 'none');
      set(hs, 'FaceLighting', 'none');

      % if facealpha is an array with number of elements equal to the number of vertices
      if size(pos,1)==numel(facealpha)
        set(hs, 'FaceVertexAlphaData', facealpha);
        set(hs, 'FaceAlpha', 'interp');
      elseif ~isempty(pos) && numel(facealpha)==1 && facealpha~=1
        % the default is 1, so that does not have to be set
        set(hs, 'FaceAlpha', facealpha);
      end

    otherwise
      ft_error('unsupported topostyle');
  end % switch contourstyle
end % plot the interpolated topography


axis off
axis vis3d
axis equal

if istrue(axes_)
  % plot the 3D axes, this depends on the units and coordsys
  ft_plot_axes([], 'coordsys', coordsys, 'unit', unit);
end

if ~isempty(coordsys)
  % add a context sensitive menu to change the 3d viewpoint to top|bottom|left|right|front|back
  menu_viewpoint(gca, coordsys)
end

if ~holdflag
  hold off
end
end

function [tri] = projecttri(pos, method)

% PROJECTTRI makes a closed triangulation of a list of vertices by
% projecting them onto a unit sphere and subsequently by constructing
% a convex hull triangulation.
%
% Use as
%   tri = projecttri(pos, method)
% where method is either 'convhull' (default) or 'delaunay'.
%
% See also SURFACE_NORMALS, PCNORMALS, ELPROJ

% Copyright (C) 2006-2019, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

tmp = pos;
tmp(:,1) = tmp(:,1) - mean(tmp(:,1));
tmp(:,2) = tmp(:,2) - mean(tmp(:,2));
tmp(:,3) = tmp(:,3) - mean(tmp(:,3));
r = rank(tmp);
switch r
    case 0
        ft_warning('vertices are lying on a single point, cannot make triangulation');
        tri = zeros(0,3);
        return
    case 1
        ft_warning('vertices are lying on a straight line, cannot make triangulation');
        tri = zeros(0,3);
        return
    case 2
        if nargin<2
            method = 'delaunay';
        end
    case 3
        if nargin<2
            method = 'convhull';
        end
    otherwise
        ft_error('unexpected input');
end

switch method
    case 'convhull'
        ori = (min(pos) + max(pos))./2;
        pos(:,1) = pos(:,1) - ori(1);
        pos(:,2) = pos(:,2) - ori(2);
        pos(:,3) = pos(:,3) - ori(3);
        nrm = sqrt(sum(pos.^2, 2));
        pos(:,1) = pos(:,1)./nrm;
        pos(:,2) = pos(:,2)./nrm;
        pos(:,3) = pos(:,3)./nrm;
        tri = convhulln(pos);
        if strcmp(surface_orientation(pos, tri), 'inward')
            % make the surface outward oriented
            tri = fliplr(tri);
        end

    case 'delaunay'
        if all(pos(:,1)==0)
            % this can happen with simulated electrode grids
            prj = pos(:,[2 3]);
        elseif all(pos(:,2)==0)
            % this can happen with simulated electrode grids
            prj = pos(:,[1 3]);
        elseif all(pos(:,3)==0)
            % this can happen with simulated electrode grids
            prj = pos(:,[1 2]);
        else
            % make a 2D triangulation of the projected points using delaunay
            prj = elproj(pos);
        end
        tri = delaunay(prj(:,1), prj(:,2));

    otherwise
        ft_error('unsupported method');
end
end

function [proj] = elproj(pos, method)

% ELPROJ makes a azimuthal projection of a 3D electrode cloud
%  on a plane tangent to the sphere fitted through the electrodes
%  the projection is along the z-axis
%
%  [proj] = elproj([x, y, z], 'method');
%
% Method should be one of these:
%     'gnomic'
%     'stereographic'
%     'orthographic'
%     'inverse'
%     'polar'
%
% Imagine a plane being placed against (tangent to) a globe. If
% a light source inside the globe projects the graticule onto
% the plane the result would be a planar, or azimuthal, map
% projection. If the imaginary light is inside the globe a Gnomonic
% projection results, if the light is antipodal a Sterographic,
% and if at infinity, an Orthographic.
%
% The default projection is a polar projection (BESA like).
% An inverse projection is the opposite of the default polar projection.
%
% See also PROJECTTRI

% Copyright (C) 2000-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

x = pos(:,1);
y = pos(:,2);
if size(pos, 2)==3
    z = pos(:,3);
end

if nargin<2
    method='polar';
end

if strcmp(method, 'orthographic')
    % this method compresses the lowest electrodes very much
    % electrodes on the bottom half of the sphere are folded inwards
    xp = x;
    yp = y;
    num = length(find(z<0));
    str = sprintf('%d electrodes may be folded inwards in orthographic projection\n', num);
    if num
        ft_warning(str);
    end
    proj = [xp yp];

elseif strcmp(method, 'gnomic')
    % the lightsource is in the middle of the sphere
    % electrodes on the equator are projected at infinity
    % electrodes below the equator are not projected at all
    rad = mean(sqrt(x.^2 + y.^2 + z.^2));
    phi = cart2pol(x, y);
    th  = atan(sqrt(x.^2 + y.^2) ./ z);
    xp  = cos(phi) .* tan(th) .* rad;
    yp  = sin(phi) .* tan(th) .* rad;
    num = length(find(th==pi/2 | z<0));
    str = sprintf('removing %d electrodes from gnomic projection\n', num);
    if num
        ft_warning(str);
    end
    xp(find(th==pi/2 | z<0)) = NaN;
    yp(find(th==pi/2 | z<0)) = NaN;
    proj = [xp yp];

elseif strcmp(method, 'stereographic')
    % the lightsource is antipodal (on the south-pole)
    rad = mean(sqrt(x.^2 + y.^2 + z.^2));
    z   = z + rad;
    phi = cart2pol(x, y);
    th  = atan(sqrt(x.^2 + y.^2) ./ z);
    xp  = cos(phi) .* tan(th) .* rad * 2;
    yp  = sin(phi) .* tan(th) .* rad * 2;
    num = length(find(th==pi/2 | z<0));
    str = sprintf('removing %d electrodes from stereographic projection\n', num);
    if num
        ft_warning(str);
    end
    xp(find(th==pi/2 | z<0)) = NaN;
    yp(find(th==pi/2 | z<0)) = NaN;
    proj = [xp yp];

elseif strcmp(method, 'inverse')
    % compute the inverse projection of the default angular projection
    [th, r] = cart2pol(x, y);
    [xi, yi, zi] = sph2cart(th, pi/2 - r, 1);
    proj = [xi, yi, zi];

elseif strcmp(method, 'polar')
    % use default angular projection
    [az, el, r] = cart2sph(x, y, z);
    [x, y] = pol2cart(az, pi/2 - el);
    proj = [x, y];

else
    ft_error('unsupported method "%s"', method);
end
end
