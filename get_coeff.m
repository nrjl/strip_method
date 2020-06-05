function [polars, Re, def] = get_coeff(filename, n_alfa, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		get_coeff
%
% PURPOSE:		Extract aerofoil coefficient data from an Excel file
%
% SYNTAX:		[polars, Re] = get_coeff(filename)
%				[polars, Re] = get_coeff(filename, n_alfa, sheets)
%
% INPUTS:		filename	- Filename and location as a string
%				n_alfa		- number of angles of attack
%				sheets		- Cell array of sheet names (surface
%								deflections)
%
% OUTPUTS:		polars		- Coefficient array (dimensions are alfa,
%								coeff, Re, def)
%				Re			- Array of Reynolds numbers
%				def			- Array of control surface deflections
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2007
%
% MODIFIED:     July 2007
%
% See also:		
%--------------------------------------------------------------------------

switch nargin
	case 1
		error('wing_model:get_coeff:n_argin', ...
			'You must enter at least the filename and number of alphas')
	case 2
		sheets = {'Sheet1'};
	case 3
		sheets = varargin{1};
	otherwise
		error('wing_model:get_coeff:n_argin', ...
			'Too many input arguments')
end

nsheets = length(sheets);

if nsheets <= 1
	def = [];
else
	def = str2double(sheets);
end

% Get Reynolds Numbers from first sheet
Re = zeros(1,20);

for j = 1:20
	loc = ['A', num2str(7 + (n_alfa+4)*(j-1))];
	[temp, Re_str] = xlsread(filename, sheets{1}, loc);
	
	if size(Re_str)
		Re_str = Re_str{1};
		Re(j) = str2double(Re_str(6:end));
	else
		break
	end
end

Re = Re(1:j-1);
nRe = j-1;

polars = zeros(n_alfa, 9, length(Re), nsheets);

for i = 1:nsheets
	for j = 1:nRe
		loc_start = ['A', num2str(10 + (n_alfa+4)*(j-1))];
		loc_stop  = ['I', num2str(5 + j*(n_alfa+4))];
		locs = [loc_start, ':', loc_stop];
		polars(:,:,j,i) = xlsread(filename, sheets{i}, locs);
	end	
end

