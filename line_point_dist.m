function d = line_point_dist(x, y, x_star, sqdist)
%--------------------------------------------------------------------------
%
% FUNCTION:		line_point_dist
%
% PURPOSE:		Find the perpendicular distance from a line to a point, 
%				with the line specified by a point (x) and a vector (y), 
%				the target point/s x_star
%               
% SYNTAX:		d = line_point_dist(x, y, x_star)
%
% INPUTS:		x	- line origin(s) - [n×d]
%				y	- line vector directions - [n×d]
%				x_star	- points - [k×d]
%
% OUTPUTS:		d	- distances [n×k]
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		March 2010
%
% MODIFIED:     March 2009
%
% See also:		
%--------------------------------------------------------------------------

x_full = permute(repmat(x, [1, 1, size(x_star, 1)]), [1, 3, 2]);


% Issue if y = [0 0 0], so that line is undefined. Maybe just use distance
% between points in this case? Somewhat computationally expensive, but fine
% for the moment
y_dist = sqrt(sum(y.^2, 2));

if nargin < 4
	sqdist = square_dist(x_star, x);
end
y_dist = y_dist + ~y_dist;		% Hack to prevent non-zero
y = y./y_dist;

y_full = permute(repmat(y, [1, 1, size(x_star, 1)]), [1, 3, 2]);
x_star_full = permute(repmat(x_star, [1, 1, size(x, 1)]), [3, 1, 2]);

d = abs(cross((x_star_full - x_full), (x_star_full - x_full - y_full))) + repmat(~y_dist, [1, size(x_star, 1)]).*sqdist;