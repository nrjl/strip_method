function uo = find_uo(chords)
%--------------------------------------------------------------------------
%
% FUNCTION:		find_uo
%
% PURPOSE:		Scaled vector length to locate zero spanwise moment point
% 				(POA = qc + uo*u) (uo is a scalar and u is vector)
%               
% INPUTS:		chords	- 2 x n matrix of chords
%
% OUTPUTS:		uo		- Scaled vector to POA
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		August 2007
%
% MODIFIED:     September 2007
%
% See also:		blade_element, strip_forces, segmenter
%--------------------------------------------------------------------------
ca = chords(1:end-1);
cb = chords(2:end);

if cb == ca
	uo = 0.5;
else
	uo = (sqrt((ca.^2 + cb.^2)/2)-ca)./(cb-ca);
end