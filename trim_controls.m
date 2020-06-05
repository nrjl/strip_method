function [X, U, varargout] = trim_controls(V, h, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		trim_controls
%
% PURPOSE:		trim the elevator control for constant airspeed flight at a
%				specified airspeed and altitude
%
% SYNTAX:		[X, U] = trim_controls(V, h)
%				[X, U] = trim_controls(V, h, X0, U0)
%				[X, U, success] = trim_controls(V, h, X0, U0)
%
% INPUTS:		V		- Target airspeed (m/s)
%				h		- Altitude (m)
%				X0		- Initial state estimate (optional)
%				U0		- Initial control estimate (optional)
%
%
% OUTPUTS:		X	- State vector
%				U	- Required control vector
%				success - Flag indicating the success of the 
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		state_rates, cruise_vel
%--------------------------------------------------------------------------
windf = 0;		% wind specified flag
w_wind = 0;		% Vertical wind component (for testing)
max_da = 2*pi/180;
max_dg = 2*pi/180;
max_de = 0.1;

if nargin == 3
	w_wind = varargin{1}; windf = 1;
elseif nargin == 4
	X = varargin{1};
	U = varargin{2};
elseif nargin == 5
	X = varargin{1};
	U = varargin{2};
	w_wind = varargin{3}; windf = 1;
elseif nargin ~= 2
	error('wing_model:trim_controls:nargin', 'Incorrect number of inputs');
end

fprintf(1, ['\nCalculating elevator trim position for %0.4g m/s TAS ', ...
	'at %0.3g m altitude\n'], V, h);

% Convergence tolerance
tol = 1e-10;

% Guess initial values
if nargin == 6
	Vw_e = wind_field(X(10:12), windf*w_wind);
	body_v(1:3) = X(1:3) - calc_Ceb(X(7:9))'*Vw_e;
	alfa = atan2(body_v(3), body_v(1));
	gamma = X(8) - alfa;
	de = U(1);
	Xbar = [alfa; gamma; de];
else
	% Guess initial values
	alfa = 1/180*pi;
	de = 0;
	gamma = -1*pi/180;	
	Xbar = [alfa; gamma; de];
	
	% Set state vector
	X = zeros(12,1);
	X(8) = alfa+gamma;
	X(10:12) = [0; 0; -h];
	U = zeros(4,1);
	U(1) = -0.1;
end

% Get wind values
Vw_e = wind_field(X(10:12));%Vw_e = wind_field(X(10:12), windf*w_wind);
X(1:3) = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;
X(12)= -h;

%% Perturbation
d_alfa = 1e-4;
d_gamma = 1e-4;
d_de = 0.01;
max_dX = [max_da; max_dg; max_de];

count = 2;
n_max = 50;
err = zeros(1,n_max);
err(1) = 2*tol;

conv = zeros(3, n_max);

while (count <= n_max) && (err(count-1) > tol)
	
	if count > 2
		X(8) = Xbar(1) + Xbar(2);
		X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + calc_Ceb(X(7:9))'*Vw_e;
		U(1) = Xbar(3);
	end
	Xdot = state_rates(X, U);
	
	% Alfa perturbation
	Xbar(1) = Xbar(1)+d_alfa;
    X(8) = Xbar(1) + Xbar(2);
	X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + calc_Ceb(X(7:9))'*Vw_e;
	U(1) = Xbar(3);
    Xdotnew = state_rates(X, U);
	J(:,1) = (Xdotnew([1,3,5]) - Xdot([1,3,5]))./d_alfa;
    Xbar(1) = Xbar(1)-d_alfa;

	% Gamma perturbation
    Xbar(2) = Xbar(2)+d_gamma;
	X(8) = Xbar(1) + Xbar(2);
	X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + calc_Ceb(X(7:9))'*Vw_e;
	U(1) = Xbar(3);
    Xdotnew = state_rates(X, U);
	J(:,2) = (Xdotnew([1,3,5]) - Xdot([1,3,5]))./d_gamma;
    Xbar(2) = Xbar(2)-d_gamma;

	% Elevator perturbation
    Xbar(3) = Xbar(3)+d_de;
	X(8) = Xbar(1) + Xbar(2);
	X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + calc_Ceb(X(7:9))'*Vw_e;
	U(1) = Xbar(3);
    Xdotnew = state_rates(X, U);
	J(:,3) = (Xdotnew([1,3,5]) - Xdot([1,3,5]))./d_de;
    Xbar(3) = Xbar(3)-d_de;

	dXbar = - inv(J)*Xdot([1,3,5]);
	Xbar2 = Xbar + dXbar;
	err(count) = (Xbar2 - Xbar)'*(Xbar2 - Xbar);
	
	% Saturate changes (to prevent unreasonable estimates)
	sat = (abs(dXbar) > max_dX);
	Xbar2 = Xbar + dXbar.*(~sat) + sign(dXbar).*max_dX.*sat;
	
	d_alfa = (Xbar2(1) - Xbar(1))/10;
	d_gamma = (Xbar2(2) - Xbar(2))/10;
	d_de = (Xbar2(3) - Xbar(3))/10;
			
	Xbar = Xbar2;
	
	conv(:,count-1) = Xbar;
	
	count = count+1;
end

Ceb = calc_Ceb(X(7:9));
X(8) = Xbar(1) + Xbar(2);
X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + Ceb'*Vw_e;
U(1) = Xbar(3);

Ve = Ceb*X(1:3);
gam = -atan2(Ve(3), Ve(1))*180/pi;

success = 1;
if count >= n_max
	fprintf(1, '\nSolution did not converge after %i iterations\n', count-2)
	figure(10);
	plot(err(1:count-1));
	success = 0;
else
	fprintf(1, '\nSolution converged after %i iterations\n', count-2)
% 	def = plane_aero(4).c_lim * [Xbar(3)> 0; Xbar(3)<0]*Xbar(3);
	fprintf(1, ['Climb angle %0.3g deg\nAngle of attack %0.3g deg\n', ...
		'Elevator = %0.3g\n'], gam, Xbar(1)*180/pi, ...
		Xbar(3));
end

if nargin == 3
	varargout{1} = success;
end

%% SCRAP
% sum(abs(Xdot([1,3]))); %sqrt(Xdot([1,3])'*Xdot([1,3]));