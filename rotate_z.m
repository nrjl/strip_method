function C = rotate_z(psi, units)
% Function to create direction cosine matrix for rotation about the z-axis

if nargin == 2
    if strcmp(units, 'deg')
        psi = psi*pi/180;
    else
        psi = -psi;
    end
end

psi = -psi;

C = [ cos(psi), sin(psi), 0;
     -sin(psi), cos(psi), 0;
         0    ,   0     , 1];