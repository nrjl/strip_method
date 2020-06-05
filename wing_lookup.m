function [Cl, Cd, Cm] = wing_lookup(alfa, Re, def, coeffs)
%--------------------------------------------------------------------------
%
% FUNCTION:		wing_lookup
%
% PURPOSE:		Look up wing coefficients based on angle of attack and 
%				known lookup data
%               
% INPUTS:		alfa	- angle of attack (deg)
%				Re		- Reynolds number
%				def		- Control surface deflection (scaled to [-1, 1])
%				coeffs	- section coefficients lookup table
%
% OUTPUTS:		Cl		- Lift coefficent
%				Cd		- Drag coefficent
%				Cm		- Pitching moment coefficent
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		August 2007
%
% MODIFIED:		September 2007
%
% See also:		strip_method, strip_forces, plot_polars
%--------------------------------------------------------------------------

%Columns:   1       2   3   4   5   6   7   8   9
%           alpha	Cl	Cd	Cm 	TU	TL	SU	SL	L/D
global PLANE_AERO

Re_range = PLANE_AERO(coeffs).profile{2};
def_range = PLANE_AERO(coeffs).profile{3};
alfa_range = transpose(PLANE_AERO(coeffs).profile{1}(:,1,1,1));

if (Re <= Re_range(1))
%  	disp('Reynolds number lower limit used')
	Re = Re_range(1);
elseif (Re >= Re_range(end))
% 	disp('Reynolds number higher limit used')
	Re = Re_range(end);
end

if (alfa < -12) || (alfa > 12)
% 	warning('blade_element:wing_lookup:alpha', ['Angle of attack %0.4g',...
% 		'is outside lookup range, using flat plate model'], alfa);
	Cmo = interp3(Re_range, alfa_range, def_range, squeeze(PLANE_AERO(coeffs).profile{1}(:,4,:,:)), Re, 10, def);
	[Cl, Cd, Cm] = flat_plate(alfa, def, 1.4, Cmo);
	return	
end

Cl = interp3(Re_range, alfa_range, def_range, squeeze(PLANE_AERO(coeffs).profile{1}(:,2,:,:)), Re, alfa, def);
Cd = interp3(Re_range, alfa_range, def_range, squeeze(PLANE_AERO(coeffs).profile{1}(:,3,:,:)), Re, alfa, def);
Cm = interp3(Re_range, alfa_range, def_range, squeeze(PLANE_AERO(coeffs).profile{1}(:,4,:,:)), Re, alfa, def);


% ---------------------------------------------------------------------
% SCRAPBOOK
%
% Cl = interp2(Re_range, coeff(:,1,1)', squeeze(coeff(:,2,:)), Re, alfa, '*linear');
% Cd = interp2(Re_range, coeff(:,1,1)', squeeze(coeff(:,3,:)), Re, alfa, '*linear');
% Cm = interp2(Re_range, coeff(:,1,1)', squeeze(coeff(:,4,:)), Re, alfa, '*linear');
%
% 	Re_ind = 1;
%
% 	Re_ind = length(Re_range);
% else
% 	Re_L = sum(Re_range < Re)
% 	Re_ind = (Re - Re_range(Re_L))/(Re_range(Re_L+1) - Re_range(Re_L));
%
% alfa_L = sum(coeff(:,1,1));
% alfa_ind = (alfa - coeff(alpha_L,1,1))/(coeff(Re_L+1,1,1)-coeff(Re_L,1,1));