function out = saturate(in, val, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		saturate
%
% PURPOSE:		Saturate a value to specified limits
%
% SYNTAX:		out = saturate(in, max)			Sets limits to [-max, max]
%				out = saturate(in, min, max)	Sets limits to [ min, max]
%
% INPUTS:		in	- Input value
%				max	- Maximum value
%				min	- Minimum value
%
% NOTE: Limit values must be either a single value or the same size as 'in'
%
% OUTPUTS:		out - Output value
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2008
%
% MODIFIED:
%
% See also:		log_profile, wind_field
%--------------------------------------------------------------------------

if nargin < 3
	min = -val;
	max = val;
else
	min = val;
	max = varargin{1};
end

u_lim = (in > max);
d_lim = (in < min);
out = in.*~(u_lim | d_lim) + u_lim.*max  + d_lim.*min;
