function C = rotate_y(theta, units)
% Function to create direction cosine matrix for rotation about the y-axis

if nargin < 2
    units = 'rad';
elseif strcmp(units, 'deg')
    theta = theta*pi/180;
else
    units = 'rad';
end

theta = -theta;

C = [ cos(theta),     0    , -sin(theta);
           0    ,     1    ,      0     ;
      sin(theta),     0    , cos(theta)];