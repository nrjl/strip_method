function X_new = runge_kutta4(X_old, U, dt)
%--------------------------------------------------------------------------
%
% FUNCTION:		runge_kutta4
%
% PURPOSE:		RK4 integration
%               
% SYNTAX:		[X_new] = runge_kutta4(X_old, U, dt)
%
% INPUTS:		X_old	- Old state vector
%				U		- Corresponding control vector
%				dt		- Time step
%
% OUTPUTS:		X_new	- New state vector
%
% AUTHOR:		Yao-Li Shen, Nicholas Lawrance
%
% CREATED:		September 2005
%
% MODIFIED:     October 2007
%
% See also:		state_rates
%--------------------------------------------------------------------------
F1 = state_rates(X_old, U);
F2 = state_rates(X_old+dt/2*F1, U);
F3 = state_rates(X_old+dt/2*F2, U);
F4 = state_rates(X_old+dt*F3, U);

X_dot = 1/6 * (F1 + 2*F2 + 2*F3 + F4);
X_new = X_old + dt*X_dot;

% X_new = X_old + dt*state_rates(X_old, U);

