function U_out = servo_check(U_in, rate, dt, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		servo_check
%
% PURPOSE:		Check that servo requests are achievable (saturation and
%				servo rate limits)
%
% SYNTAX:		U_out = servo_check(U_in, rate)
%				U_out = servo_check(U_in, rate, lim)
%
% INPUTS:		U_in	- [n×2] matrix of n control inputs (previous and
%							requested timestep in respective columns)
%				rate	- servo rate limits (measured in units per second)
%				dt		- time step
%				lim		- servo limit (default 1)
%
% OUTPUTS:		(display to command window)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		
%--------------------------------------------------------------------------

% Check if limits are specified
if nargin == 4
	lim = varargin{1};
else
	lim = 1;
end

% Check rate limits
U_dot = (U_in(:,2) - U_in(:,1))/dt;
U_break = abs(U_dot) > (rate*ones(size(U_dot)));
U_lim = U_in(:,1) + sign(U_dot)*rate*dt;
U_out = ~U_break.*U_in(:,2) + U_break.*U_lim;


% Saturate limits
U_u = (U_out > lim);
U_d = (U_out < -lim);
U_out = U_out.*~(U_u | U_d) + U_u  - U_d;

