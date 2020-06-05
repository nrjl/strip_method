function [poa_e, F, M, varargout] = strip_forces(qc, x, ai, wing_handle, coeffs, def, mirror)
%--------------------------------------------------------------------------
%
% FUNCTION:		strip_forces
%
% PURPOSE:		Calculate strip forces on a single surface
%               
% SYNTAX:		[poa_e, F, M] = strip_forces(qc, x, ai, wing_handle, coeff,
%					Re_coeff, def)
%
% INPUTS:		qc	- vector of quarter chord points in body axes (m)
%				x	- state vector
%				ai	- base angle of incidence (deg)
%				coeffs		- lookup table of wing coefficients
%				def			- control surface deflection (deg)
%				mirror		- flag if surface is mirrored (1 for mirrored,
%								else 0)
%
% OUTPUTS:		poa_e	- Points of force application in earth frame
%				F		- 2D force vector in wing coordinates
%				M		- Total moment around CG provided by each section
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		August 2007
%
% MODIFIED:		September 2007
%
% See also:		strip_method, segmenter, wind_field
%--------------------------------------------------------------------------

if ~isa(wing_handle, 'function_handle')
    error('blade_element:strip_forces:nothandle',...
        ['Fourth input argument must be a function handle \n',...
        'for the aero coefficients as a function of incidence'])
end

n_poa = size(qc, 2)-1;			% Number of points of application

F = zeros(3, n_poa);			% Initialise forces and moments
M = zeros(3, n_poa);
stalled = zeros(1, n_poa);

[T, P, rho, a, mu, nu] = atmos(-x(12), 0);		% Atmospheric properties

u = diff(qc(1:3,:), 1, 2);						% Vectors between each qc point
cbar = (qc(4,1:end-1) + qc(4,2:end))/2;			% Mean chord (assume linear taper)
d = sqrt(u(1,:).^2 + u(2,:).^2 + u(3,:).^2);	% Strip widths
uo = find_uo(qc(4,:));							% Get vector scales to zero roll POA
poa = qc(1:3,1:end-1) + ([1;1;1]*uo).*u;		% Points of application

Ceb = calc_Ceb(x(7:9));	% Earth to body DCM
Cbe = Ceb';

V_cgg_b = x(1:3);	% Velocity of CG with respect to ground in body axes
omega_b = x(4:6);	% Rotation rates of body

poa_e = Ceb*poa + x(10:12)*ones(1, size(poa, 2));
V_wg_e = wind_field(poa_e);	% Velocity of wind with respect to ground in earth axes
V_wg_b = Cbe*V_wg_e;		% Velocity of wind with respect to ground in body axes

% lin_flag = 0;				% Flag to note that local linear alpha approximation is being used

% Finite wing vortex correction
b = sqrt(sum((qc(2:3,end) - qc(2:3,1)).^2))*(1 + mirror);
AR = b*b/((cbar*d')*(1+mirror));
oswald_eff = 0.9;

% Forces loop for each strip (applied at points of application)
for i = 1:n_poa
    V_pw_b = V_cgg_b + cross(omega_b, poa(:,i)) - V_wg_b(:,i); % Velocity of p wrt wind in body
	if u(2,i)
		lr = sign(u(2,i));						% Left or right surface
	else
		lr = 1;
	end
	
	dihedral = lr*atan2(u(3,i),abs(u(2,i)));	% Effective dihedral of section
    
	u_cone = rotate_x(-dihedral)*u(:,i);
	sweep = lr*atan2(u_cone(1),abs(u_cone(2)));	% Effective sweep of section
    	
% 	if abs(abs(dihedral)- pi/2) < 1e-7; sweep = 0; end;
	
	% Body to wing DCM
	Cwb = rotate_y(-ai(i)*pi/180)*rotate_z(sweep)*rotate_x(-dihedral);
	V_pw_w = Cwb*V_pw_b;	% Velocity of point wrt wind in wing axes
    
    alfa_i = atan2(V_pw_w(3), V_pw_w(1));		% Effective incidence (rad)
	if alfa_i > 7*pi/180 || alfa_i < -10*pi/180
		stalled(i) = 1;
	end
	
    qbar = 0.5*rho*(V_pw_w(3)^2 + V_pw_w(1)^2);	% Dynamic pressure
	Re_i = sqrt(V_pw_w(3)^2 + V_pw_w(1)^2)*qc(4,i)/nu;	% Reynolds
    
	% This section uses a linear approximation for the coefficients if the
	% angle of attack is not near stall
% 	if i == 1 && (abs(alfa_i < 6*pi/180) || abs(alfa_i > 12*pi/180))
		[Cl, Cd, Cm] = wing_handle(alfa_i*180/pi, Re_i, def(i), coeffs);	% Coefficients at alpha
% 		[Cl2, Cd2, Cm2] = wing_handle(alfa_i*180/pi+0.1, Re_i, def(i), coeffs);	% Coefficients at alpha + 0.5
% 		al0 = alfa_i + 0.05*pi/180;
% 		Cl0 = (Cl2 + Cl)/2; dcl_da = (Cl2 - Cl)/(0.1*pi/180);
% 		Cd0 = (Cd2 + Cd)/2; dcd_da = (Cd2 - Cd)/(0.1*pi/180);
% 		Cm0 = (Cm2 + Cm)/2; dcm_da = (Cm2 - Cm)/(0.1*pi/180);
% 		lin_flag = 1;
% 	else
% 		if lin_flag
% 			Cl = Cl0 + dcl_da*(alfa_i - al0);
% 			Cd = Cd0 + dcd_da*(alfa_i - al0);
% 			Cm = Cm0 + dcm_da*(alfa_i - al0);
% 		else
% 			[Cl, Cd, Cm] = wing_handle(alfa_i*180/pi, Re_i, def(i), coeffs);
% 		end
% 	end
	Cd = Cd + Cl*Cl/(pi*AR*oswald_eff);
	
    Li = Cl*qbar*cbar(i)*d(i);					% Section lift
    Di = Cd*qbar*cbar(i)*d(i);					% Section drag
    Mi = Cm*qbar*cbar(i)^2*d(i);				% Section pitching moment
    
    % Forces in wing axes
    Fi_w = [Li*sin(alfa_i) - Di*cos(alfa_i);...
            0;...
            -Li*cos(alfa_i) - Di*sin(alfa_i)];
    
    % Restore forces and moments to body axes
    F(:,i) = Cwb'*Fi_w;
    M(:,i) = Cwb'*[0;Mi;0] + cross(poa(:,i), F(:,i));
end

if nargout == 4
	varargout(1) = {stalled};
end