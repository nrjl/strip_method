function [Cl, Cd, Cm] = flat_plate(alfa, def, Cl_max, Cm_o)
%--------------------------------------------------------------------------
%
% FUNCTION:		flat_plate
%
% PURPOSE:		calcualte modified aerodynamic coefficients of a flat plate
%               
% SYNTAX:		[Cl, Cd, Cm] = flat_plate(alfa, def, Cl_max, Cm_o)
%
% INPUTS:		alfa	- Angle of attack (degrees)
%				def		- Control surface deflection (deg)
%				Cl_max	- Maximum Cl for the surface
%				Cm_o	- Pitch moment coefficient at zero angle of attack
%
% OUTPUTS:		Cl		- Lift coefficient
%				Cd		- Drag coefficient
%				Cm		- Pitch moment coefficient
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2005
%
% MODIFIED:     October 2007
%
% See also:		aero_plot, strip_plot
%--------------------------------------------------------------------------

Cl = (Cl_max+sign(alfa)*def/20 - 0.15*(alfa<0))*sind(2*alfa);

Cd_max = 2*(Cl_max + sign(alfa)*def/20);
Cd = Cd_max*(sind(alfa).^2);

Cm = -Cd/4.*sind(alfa) + Cm_o;