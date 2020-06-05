function C = rotate_x(phi, units)
% Function to create direction cosine matrix for rotation about the x-axis

if nargin < 2
    units = 'rad';
elseif strcmp(units, 'deg')
    phi = phi*pi/180;
else
    units = 'rad';
end

phi = -phi;

C = [ 1,     0    ,   0     ;
      0,  cos(phi), sin(phi);
      0, -sin(phi), cos(phi)];