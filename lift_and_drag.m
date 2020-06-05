function [L, D, CL, CD, Cd0, n] = lift_and_drag(varargin)
% Lift and drag calculator for SBXC

global g
global PLANE_AERO PLANE_PARAM
e = 0.85;
	
if nargin == 0
	V = 15;
	alt = 100;


	if isempty(g); initialise_SBXC; end
	[X0, U0] = trim_controls(V, alt);
else
	X0 = varargin{1};
	U0 = varargin{2};
end

%%
[S, AR] = plane_properties(PLANE_AERO, PLANE_PARAM);
V = sqrt(sum(X0(1:3).^2));
[T, P, rho] = atmos(-X0(12), 0);
[FA, MA] = strip_method(X0, U0);
[FF, MF] = fuselage_forces(X0);
F = FA + FF;

alfa = atan2(X0(3), X0(1));

L = -F(3)*cos(alfa) + F(1)*sin(alfa);
D = -F(3)*sin(alfa) - F(1)*cos(alfa);

CL = L/(0.5*rho*V*V*S);
CD = D/(0.5*rho*V*V*S);

Cd0 = CD - CL*CL/(pi*AR*e); %- polyval(p_Cd, CL)

n = L/(PLANE_PARAM.m*g);

fprintf('\n  L =%9.5g\t  D =%9.5g\n CL =%9.5g\t CD =%9.5g\nCd0 =%9.5g\n',...
	L, D, CL, CD, Cd0);
fprintf('Wing loading = %0.5g\n\n', n);


