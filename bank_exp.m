function bank_out = bank_exp(K, V_min, V_in)
%--------------------------------------------------------------------------
%
% FUNCTION:		bank_exp
%
% PURPOSE:		Calculate the minimum energy bank angle for a turn at a
%				given airspeed
%
% SYNTAX:		bank_out = bank_exp(K, V_min, V_in)
%
% INPUTS:		K		- Exponential function values (obtained using
%							Matlab's non-linear optimisation toolkit)
%				V_min	- Minimum allowable turn speed
%				V_in	- Current airspeed
%
% OUTPUTS:		bank_out- Energy optimum bank angle in rads
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		May 2008
%
% MODIFIED:
%
% See also:		
%--------------------------------------------------------------------------
bank_out = K(1)*(1 - K(2).^(V_min - V_in));
bank_out = bank_out.*(bank_out > 0)*pi/180;