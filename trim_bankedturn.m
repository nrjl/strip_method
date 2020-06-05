function [X, U, varargout] = ...
	trim_bankedturn(V, bank, h, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		trim_bankedturn
%
% PURPOSE:		trim aileron, elevator and rudder for a zero-sideslip
%					steady banked turn
%
% SYNTAX:		[X, U] = trim_controls(V, bank, h)
%				[X, U, fail_flag] = trim_controls(V, bank, h, X, U)
%
% INPUTS:		V	- Target airspeed (m/s)
%				bank- Bank attitude (deg)
%				h	- Altitude (m) (default value is 300m if not specified)
%
% OUTPUTS:		X	- State vector
%				U	- Required control vector
%				fail_flag - Flag determining success of operation (1 is
%								fail)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		state_rates, cruise_vel, trim_controls
%--------------------------------------------------------------------------
global PLANE_PARAM

if nargin == 5
	X_guess = varargin{1};
	U_guess = varargin{2};
elseif nargin ~= 3
	error('wing_model:trim_controls:nargin', 'Incorrect number of inputs');
end

fprintf(1, ['\nCalculating trim controls for coordinated turn at ', ...
	'%0.4g m/s TAS, %0.3g deg bank, at %0.3g m altitude\n'], V, bank, h);

% Convergence tolerance
tol = 1e-10*(pi/180)^2;
bank = bank*pi/180;


% Check if initial estimate has been provided
if nargin == 7
	Vw_e = wind_field(X_guess(10:12));
	body_v(1:3) = X_guess(1:3) - calc_Ceb(X_guess(7:9))'*Vw_e;
	alfa = atan2(body_v(3),body_v(1));
	Xbar = [U_guess(1:3); X_guess(6); X_guess(8); alfa];
	X = X_guess;
else
	% Guess initial values
	de = -0.0863;
	da = 0.2560;
	dr = 0.1126;
	r = 0.3040;
	theta = -0.0034;
	alfa = 1.1325*pi/180;
	Xbar = [de; da; dr; r; theta; alfa];

	% Set state vector
	X = zeros(12,1);
	X(4) = -r*tan(theta)*(sin(bank)*tan(bank) + cos(theta));
	X(5) = r*tan(bank);
	X(6) = r;
	X(8) = theta;
end

X(12)= -h;
X(7) = bank;	

Vw_e = wind_field(X(10:12));
X(1:3) = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;

%% Perturbation
perturb = 1e-4*ones(6,1);
kick = 0;

count = 2;
nmax = 30;
err = zeros(1,nmax);
err(1) = 2*tol;

conv = zeros(length(Xbar), nmax);

while (count <= nmax) && (err(count-1) > tol)
	U_c = Xbar(1:3);
	U_u = (U_c >  1-perturb(1:3));
	U_d = (U_c < -1-perturb(1:3));
	Xbar(1:3) = U_c.*~(U_u | U_d) + U_u.*(U_u-perturb(1:3)) + U_d.*(-U_d-perturb(1:3));
	
	de		= Xbar(1);
	da		= Xbar(2);
	dr		= Xbar(3);
	r		= Xbar(4);
	theta	= Xbar(5);
	alfa	= Xbar(6);
	
	if count > 2
		X(1:3)   = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;
		X(4) = -r*tan(theta)*(sin(bank)*tan(bank) + cos(theta));
		X(5) = r*tan(bank);
		X(6) = r;
		X(7) = bank;
		X(8) = theta;
	end
	
	Xold = X;
	U = [de; da; dr; 0];	
	Xdot = state_rates(X, U);
	
	% Elevator
	U(1) = de + perturb(1);
	Xdotnew = state_rates(X, U);
	J(:,1) = (Xdotnew(1:6) - Xdot(1:6))./perturb(1);
    U(1) = de;
	
	% Aileron
	U(2) = da + perturb(2);
	Xdotnew = state_rates(X, U);
	J(:,2) = (Xdotnew(1:6) - Xdot(1:6))./perturb(2);
	U(2) = da;
	
	% Rudder
	U(3) = dr + perturb(3);
	Xdotnew = state_rates(X, U);
	J(:,3) = (Xdotnew(1:6) - Xdot(1:6))./perturb(3);
	U(3) = dr;
	
	% r
    r = r + perturb(4);
	X(4) = -r*tan(theta)*(sin(bank)*tan(bank) + cos(theta));
	X(5) = r*tan(bank);
	X(6) = r;
	Xdotnew = state_rates(X, U);
	J(:,4) = (Xdotnew(1:6) - Xdot(1:6))./perturb(4);
	r = Xbar(4);
	X = Xold;
	
	% theta
	theta = theta + perturb(5);
	X(4) = -r*tan(theta)*(sin(bank)*tan(bank) + cos(theta));
	X(8) = theta;
	Xdotnew = state_rates(X, U);
	J(:,5) = (Xdotnew(1:6) - Xdot(1:6))./perturb(5);
	X = Xold;
	
	% alpha
	alfa = alfa + perturb(6);
	X(1:3)   = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;
	Xdotnew = state_rates(X, U);
	J(:,6) = (Xdotnew(1:6) - Xdot(1:6))./perturb(6);
	X = Xold;

	Xbar2 = Xbar - inv(J)*Xdot(1:6);
	err(count) = (Xbar2 - Xbar)'*(Xbar2 - Xbar);	
	
	perturb = (Xbar2 - Xbar)/10;
	
	Xbar = Xbar2;	
	
	if (abs(mod(count-1, 10)) < 1e-8) || ((err(count) > err(count-1)) && count > 5) && kick == 0
		perturb = perturb.*(randn(size(perturb))/2 + 1);
		Xbar = Xbar.*(randn(size(Xbar))/5 + 1);
		kick = 1;
		% Xbar = mean(conv(:,count-20:count-2),2);
	elseif kick == 4
		kick = 0;
	else
		kick = kick+1;
	end
	
	conv(:,count-1) = Xbar;
	
	count = count+1;
end

if count >= nmax
	fprintf(1, '\nSolution did not converge after %i iterations\n',count-2)
	figure(10);
	plot(err(1:count-1));
	fail_flag = 1;
	X = NaN(12,1);
	U = NaN(4,1);
	
	if nargout == 4
		varargout(1) = {fail_flag};
		varargout(2) = {zeros(size(X))};
	end
else
	fprintf(1, '\nSolution converged after %i iterations\n', count-2)
% 	def = PLANE_AERO(4).c_lim * [Xbar(3)> 0; Xbar(3)<0]*Xbar(3);
	fprintf(1, ['Elevator = %0.3g\nAileron = %0.3g\nRudder = %0.3g\n',...
		'Body yaw = %0.3g deg/s\nPitch angle = %0.3g deg\n',...
		'Alpha = %0.3g\n'], Xbar(1), Xbar(2), Xbar(3), Xbar(4)*180/pi, ...
		Xbar(5)*180/pi, Xbar(6)*180/pi);
	fprintf(1, '\nPower loss = %0.5g W\n', Xdot(12)*9.81*PLANE_PARAM.m);
	fprintf(1, 'Turn efficiency = %0.5g deg/m\n', Xdot(9)*180/pi/Xdot(12));
	fail_flag = 0;

	r		= Xbar(4);
	theta	= Xbar(5);
	alfa	= Xbar(6);
	X(1:3)  = [V*cos(alfa); 0; V*sin(alfa)] + calc_Ceb(X(7:9))'*Vw_e;
	X(4) = -r*tan(theta)*(sin(bank)*tan(bank) + cos(theta));
	X(5) = r*tan(bank);
	X(6) = r;
	X(7) = bank;
	X(8) = theta;
	U = [Xbar(1:3); 0];
	
	if nargout == 4
		varargout(2) = {Xdot};
	end
end

if nargout >= 3
	varargout(1) = {fail_flag};
end

%% SCRAP
% sum(abs(Xdot([1,3]))); %sqrt(Xdot([1,3])'*Xdot([1,3]));