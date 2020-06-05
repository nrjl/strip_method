function xf = newton_sol(fn_handle, x0, dx, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		newton_sol
%
% PURPOSE:		Use Newton's solution to numerically solve an equation
%
% SYNTAX:		xf = newton_sol(fn_handle, x0, dx)
%				xf = newton_sol(fn_handle, x0, dx, yt)
%				xf = newton_sol(fn_handle, x0, dx, yt, dfn_handle)
%
% INPUTS:		fn_handle	- Handle to primary function
%				x0			- Initial solution estimate
%				dx			- Initial step size
%				yt			- Target function value (default zero)
%				dfn_handle	- Derivative function handle (will be solved
%					using numeric derivatives if this is not specified)
%
% OUTPUTS:		xf			- Function solution
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2007
%
% MODIFIED:     July 2007
%
% See also:		
%--------------------------------------------------------------------------
err = 1e-7;             % convergence error
yt = 0;                 % target y (default 0)
nmax = 100;             % max number of iterations

if nargin == 3
    der = 0;            % derivative flag
elseif nargin == 4
    der = 0;
    yt = varargin{1};   % target function value
elseif nargin == 5
    der = 1;
    yt = varargin{1};   % target function value
    df = varargin{2};   % derivative function handle
else
    error('nick:newton_sol:n_inputs', 'Incorrect number of inputs')
end

i = 1;
% x(1) = x0;


while (abs(dx) > err) && (i < nmax)
    y1 = fn_handle(x0);
	
    if der == 0
        y2 = fn_handle(x0+dx);
        dfdx = (y2-y1)/dx;    
    else
        dfdx = df(x0);
    end
    
    dx = (y1 - yt)/dfdx;
    
    x0 = x0 - dx;
	i = i+1;
%     x(i) = x0;   
end

if (i >= nmax)
	warning('nick:newton_sol:no_converge', 'Solution did not converge')
end

xf = x0;    