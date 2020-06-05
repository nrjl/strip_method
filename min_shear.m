function [gamma, beta, varargout] = min_shear(V, alt, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		min_shear
%
% PURPOSE:		Calculate the minimum vertical shear rate required for
%				velocity sustain flight at a given airspeed and altitiude. 
%				Climb angle can be specified. If it is not specified, the 
%				minimum required value is calculated.
%
% SYNTAX:		[gamma, beta] = min_shear(V, alt)
%				[gamma, beta] = min_shear(V, alt, gamma)
%				[gamma, beta, Edot] = min_shear(V, alt, gamma)
%
% INPUTS:		V		- Airspeed (m/s)
%				alt		- Altitude (m)
%				gamma	- Cimb angle (rad)
%
% OUTPUTS:		gamma	- Cimb angle (rad)
%				beta	- Required vertical shear rate (1/s)
%				Edot	- Energy gain rate ( = V*sin(gamma))
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		January 2008
%
% MODIFIED:     February 2008
%
% See also:		min_shear_plot, min_shear_plot3
%--------------------------------------------------------------------------
%% SETUP
global PLANE_PARAM
if nargin == 2
	find_gamma = 1;
elseif 	nargin == 3
	find_gamma = 0;
	gamma = varargin{1};
else
	error('Incorrect number of inputs');
end

m = PLANE_PARAM.m;
g = 9.81;

% Calculate body forces for specified trim condition
[X, U] = trim_controls(V, alt);
[F, M] = strip_method(X, U);

alfa = atan2(X(3), X(1));				% Calculate angle of attack
D = -F(3)*sin(alfa) - F(1)*cos(alfa);	% Calculate drag

if find_gamma
	% If ideal climb angle is required
	% Solution of a cubic equation for one real root
	
	% Solve for the climb angle, gamma
	y = (m*g/D)^(1/3);
	u = y + 2/(3*y);
	gamma = asin(1/u);
end

% Calculate required wind shear to maintain velocity
beta = (D + m*g.*sin(gamma))./(m*V.*sin(gamma).*cos(gamma));

% Calculate wind power in specified conditions
Edot_wind = m*g*V*sin(gamma);

if nargout == 3
	varargout(1) = {Edot_wind};
elseif nargout == 4
	varargout(1) = {Edot_wind};
	varargout(2) = {D};
elseif nargout > 4
	disp('Warning: too many output arguments');
end