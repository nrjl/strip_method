function X_dot = state_rates(X, U)
%--------------------------------------------------------------------------
%
% FUNCTION:		state_rates
%
% PURPOSE:		calculate aircraft state rates
%               
% SYNTAX:		Xdot = state_rates(X, U)
%
% INPUTS:		X		- state vector
%				U		- control vector
%
% OUTPUTS:		Xdot	- state rates
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:      2007
%
% See also:		strip_method, strip_forces, SBXC_def, test
%--------------------------------------------------------------------------
global g PLANE_PARAM

m = PLANE_PARAM.m;
Ixx = PLANE_PARAM.Ixx;
Iyy = PLANE_PARAM.Iyy;
Izz = PLANE_PARAM.Izz;
Ixz = PLANE_PARAM.Ixz;


X = X(:);
U = U(:);

% States
u	= X(1,1);
v	= X(2,1);
w	= X(3,1);
p	= X(4,1);
q	= X(5,1);
r	= X(6,1);
phi	= X(7,1);
the	= X(8,1);
% psi	= X(9,1);
% xe	= X(10,1);
% ye	= X(11,1);
ze	= X(12,1);

% Controls
% de = U(1,1);
% da = U(2,1);
% dr = U(3,1);
% df = U(4,1);

Ceb =  calc_Ceb(X(7:9));

% Aerodynamic forces
[AForces, AMoments] = strip_method(X, U);
[FForces, FMoments] = fuselage_forces(X);

FxA = AForces(1) + FForces(1);
FyA = AForces(2) + FForces(2);
FzA = AForces(3) + FForces(3);

% Propulsive forces
TForces = [0;0;0];
FxT = TForces(1);
FyT = TForces(2);
FzT = TForces(3);

% Moments about the body axes
L = AMoments(1) + FMoments(1);
M = AMoments(2) + FMoments(2);
N = AMoments(3) + FMoments(3);

%State Rates  - non-linear solutions
C0 = Ixx*Izz-(Ixz^2);             % required constants
C1 = Izz/C0;
C2 = Ixz/C0;
C3 = C2*(Ixx - Iyy + Izz);
C4 = C1*(Iyy-Izz)-(C2*Ixz);
C5 = 1/Iyy;
C6 = C5*Ixz;
C7 = C5*(Izz-Ixx);
C8 = Ixx/C0;
C9 = C8*(Ixx-Iyy) + (C2*Ixz);

% Gravitational force in body axes
k_ground = 10;
g_b = Ceb'*[0, 0, g*(ze<=0) + (ze>0)*-ze*k_ground]';


% Calculation of derivatives
u_dot = FxA/m + FxT/m + r*v - q*w + g_b(1);		%calculation of ud
v_dot = FyA/m + FyT/m - r*u + p*w + g_b(2);		%calculation of vd
w_dot = FzA/m + FzT/m + q*u - p*v + g_b(3);		%calculation of wd
p_dot = C3*p*q + C4*q*r + C1*L + C2*N;			%calculation of pd
q_dot = C7*p*r - C6*((p^2)-(r^2)) + C5*M;		%calculation of qd
r_dot = C9*p*q - C3*q*r + C2*L + C8*N;			%calculation of rd

phi_dot = p + q*sin(phi)*tan(the) + r*cos(phi)*tan(the);
the_dot = q*cos(phi) - r*sin(phi);
psi_dot = q*sin(phi)*sec(the) + r*cos(phi)*sec(the);


% Earth frame velocity

v_earth = Ceb*[u v w]';
xe_dot = v_earth(1,1);
ye_dot = v_earth(2,1);
ze_dot = v_earth(3,1);

X_dot = [u_dot;		v_dot;	w_dot; ...
		 p_dot;		q_dot;	r_dot; ...
		 phi_dot;	the_dot;psi_dot;...
		 xe_dot;	ye_dot; ze_dot];