function [X, U] = trim_shear(V, h, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		trim_shear
%
% PURPOSE:		trim the elevator control for constant airspeed flight in a
%				shear layer
%
% SYNTAX:		[X, U] = trim_shear(V, h)
%				[X, U] = trim_shear(V, h, X0, U0)
%
% INPUTS:		V	- Target airspeed (m/s)
%				h	- Altitude (m)
%				X0	- Initial state estimate
%				U0	- Initial control estimate
%
%
% OUTPUTS:		X	- State vector
%				U	- Required control vector
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		February 2008
%
% MODIFIED:     February 2008
%
% See also:		state_rates, cruise_vel
%--------------------------------------------------------------------------

if nargin == 4
	X = varargin{1};
	U = varargin{2};
else
	error('wing_model:trim_controls:nargin', 'Incorrect number of inputs');
end

fprintf(1, ['\nCalculating elevator trim position for %0.4g m/s TAS ', ...
	'at %0.3g m starting altitude\n'], V, h);

% Convergence tolerance
tol = 1e-10;

% Guess initial values
if nargin == 4
	Vw_e = wind_field(X(10:12));
	body_v(1:3) = X(1:3) - calc_Ceb(X(7:9))'*Vw_e;
	alfa = atan2(body_v(3), body_v(1));
	gamma = X(8) - alfa;
	de = U(1);
	Xbar = [alfa; gamma; de];
else
	% Guess initial values
	alfa = 1/180*pi;
	de = -0.1;
	gamma = 15*pi/180;	
	Xbar = [alfa; gamma; de];
	
	% Set state vector
	X = zeros(12,1);
	X(8) = alfa+gamma;
	X(10:12) = [0; 0; -h];
	U = zeros(4,1);
	U(1) = -0.1;
end

% Get wind values
Vw_e = wind_field(X(10:12));
X(1:3) = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;
X(12)= -h;

%% Perturbation
d_alfa = 1e-4;
d_gamma = 1e-4;
d_de = 0.01;

count = 2;
err = zeros(1,100);
err(1) = 2*tol;

conv = zeros(3, 100);

while (count <= 100) && (err(count-1) > tol)
	
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

	Xbar2 = Xbar - inv(J)*Xdot([1,3,5]);
	err(count) = (Xbar2 - Xbar)'*(Xbar2 - Xbar);	
	
	d_alfa = (Xbar2(1) - Xbar(1))/10;
	d_gamma = (Xbar2(2) - Xbar(2))/10;
	d_de = (Xbar2(3) - Xbar(3))/10;
			
	Xbar = Xbar2;	
	
% 	if ~mod(count-1, 20)
% 		Xbar = mean(conv(:,count-20:count-2),2);
% 	end
	
	conv(:,count-1) = Xbar;
	
	count = count+1;
end

Ceb = calc_Ceb(X(7:9));
X(8) = Xbar(1) + Xbar(2);
X(1:3) = [V*cos(Xbar(1)); 0; V*sin(Xbar(1))] + Ceb'*Vw_e;
U(1) = Xbar(3);

Ve = Ceb*X(1:3);
gam = -atan2(Ve(3), Ve(1))*180/pi;

if count >= 100
	fprintf(1, '\nSolution did not converge after %i iterations\n', count-2)
	figure(10);
	plot(err(1:count-1));
else
	fprintf(1, '\nSolution converged after %i iterations\n', count-2)
% 	def = plane_aero(4).c_lim * [Xbar(3)> 0; Xbar(3)<0]*Xbar(3);
	fprintf(1, ['Climb angle %0.3g deg\nAngle of attack %0.3g deg\n', ...
		'Elevator = %0.3g\n'], gam, Xbar(1)*180/pi, ...
		Xbar(3));
end



%% SCRAP
% sum(abs(Xdot([1,3]))); %sqrt(Xdot([1,3])'*Xdot([1,3]));