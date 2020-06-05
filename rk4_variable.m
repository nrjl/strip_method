function [X_new, lvl] = rk4_variable(X_old, U, dt, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		rk4_variable
%
% PURPOSE:		variable time step RK4 integration
%               
% SYNTAX:		X_new = rk4_variable(X_old, U, dt, tol)
%				X_new = rk4_variable(X_old, U, dt, tol, X_dt, lvl, step_check)*
%
%				*Note that only the first or second syntaxes should be used
%				from higher level functions, the X_dt specified call is for
%				recursive operation.
%
% INPUTS:		X_old	- Old state vector
%				U		- Corresponding control vector
%				dt		- Time step
%				aero	- aircraft aerodynamic definition
%				param	- aircraft parameters definition
%				tol		- convergence tolerance
%				X_dt	- single time step rk4 solution
%				lvl		- current integration level (1 is root level)
%
% OUTPUTS:		X_new	- New state vector
%				Xdot	- Single-step derivative state vector 
%							(=(X_new - X_old)/dt)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		state_rates, runge_kutta4
%--------------------------------------------------------------------------

if nargin == 7
	tol = varargin{1};
	X_dt = varargin{2};
	if isempty(X_dt); X_dt = runge_kutta4(X_old, U, dt); end;		
	lvl = varargin{3};
	step_check = varargin{4};	% 1 for upwards check, 0 for no check, -1 for down check
elseif nargin == 4
	tol = varargin{1};
	
	lvl = 1;
	step_check = -1;
else
	tol = 1e-4;
	X_dt = runge_kutta4(X_old, U, dt);
	lvl = 1;
	step_check = -1;
end

if step_check == -1			% Check if step is fine enough
	% Time halving approach
	X_new1 = runge_kutta4(X_old, U, dt/2);
	X_new2 = runge_kutta4(X_new1, U, dt/2);

	err = (X_dt - X_new2)'*(X_dt - X_new2);

	if (err > tol) && (err > 1e-6) && (lvl < 8)
		% ERROR CHECK FAILED
		[X_new, lvl] = rk4_variable(X_old, U, dt/2, tol*0.75, X_new1, lvl+1, -1);
	else
		% ERROR CHECK PASSED
		fprintf(1, 'Level %i at dt = %0.3g\n', lvl, dt/2)
		X_new = zeros(size(X_old,1), 2^lvl);
		X_new(:,1:2) = [X_new1, X_new2];

		for i = 3:2^lvl
			X_new(:,i) = runge_kutta4(X_new(:,i-1), U, dt/2);
		end		
	end
	
elseif step_check == 1 % Check if step is too fine
	X_new2 = runge_kutta4(X_dt, U, dt);
	
	X_big = runge_kutta4(X_old, U, dt*2);
	
	err = (X_big - X_new2)'*(X_big - X_new2);
	
	if ((err < tol/0.75) && (lvl > 1)) || (err < 1e-6/.75)
		% ERROR SUFFICIENTLY LOW TO MOVE UP LEVEL
		fprintf(1, 'Level %i at dt = %0.3g\n', lvl-1, dt)
% 		[X_new, lvl] = rk4_variable(X_old, U, dt*2, tol/0.75, X_dt, lvl-1, 1);

		lvl = lvl-1;
		X_new = zeros(size(X_old,1), 2^lvl);
		X_new(:,1:2) = [X_dt, X_new2];

		for i = 3:2^lvl
			X_new(:,i) = runge_kutta4(X_new(:,i-1), U, dt);
		end
	else
		% ERROR CHECK FOR LOWER LEVEL
		[X_new, lvl] = rk4_variable(X_old, U, dt, tol, X_dt, lvl, -1);
	end
	
elseif step_check == 0
	X_new = X_dt;
end

% if nargout == 2
% 	varargout(1) = {X_dot};
% end

% ----- SCRAP ----- %
% err = sum(sum(abs(diff([F1(1:9), F2(1:9), F3(1:9), F4(1:9)], 1, 2))));
% if (err > dt/100) && (dt > 0.005)
% 	fprintf(1, 'RK4 failed, attempting timestep = %0.3g\n', dt/2)
% 	X_new = runge_kutta4(X_old, U, dt/2);
% 	X_new = runge_kutta4(X_new, U, dt/2);
% end
% X_new = X_old + dt*state_rates(X_old, U);