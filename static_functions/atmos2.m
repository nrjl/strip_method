function [T, varargout] = atmos2(z, show)
%--------------------------------------------------------------------------
%
% FUNCTION:		atmos
%
% PURPOSE:		Matlab code to calculate the atmospheric properties as a
%               function of altitude. ICAO 1993 standard
%               model is used. Available for altitudes 0 - 86000m.
%               
% SYNTAX:		[T, P, rho, a, mu, nu] = atmos(z, show)
%
% INPUTS:		z - geometric altitude (m)
%               show - flag to display output (default is off)
%
% OUTPUTS:		T   - Temperature (K)
%               P   - Pressure (Pa)
%               rho - Density (kg/m^3)
%               a   - Local speed of sound (m/s)
%               mu  - Dynamic viscosity (Pa.s)
%               nu  - Kinematic viscosity (m^2/s)
%
% AUTHOR:		Nicholas Lawrance
%
% DATE:			May 2007
%
% MODIFIED:		July 2011 - ICAO 1993 Standard
%
% See also:		atmos_plot
%--------------------------------------------------------------------------

if nargin < 2
    show = 0;
end

% Constants
To = 288.15;    % K
Po = 101325;    % Pa
R = 8.31432;    % N.m/mol.K
M = 0.02896442; % kg/mol (air)
gamma = 1.40;   % (dimless)
Ra = R/M;       % J/kg.K

E = 6356766;	% Radius of Earth (m)
g0 = 9.80665;   % Gravitational acceleration at SL
GM_R = g0/Ra;	% GM/R term

% Sutherland's equation reference values for air
mu_o = 18.27e-6;% Reference viscosity
T_o  = 291.15;  % Reference temperature
C    = 120;     % Sutherland's constant

% Get Lapse rate dependence on altitude
h_layers = [0, 11, 20, 32, 47, 51, 71, 84.852]*1e3; % Layer heights in m
L_layers = [-6.5, 0, 1.0, 2.8, 0, -2.8, -2.0]/1e3;  % Layer lapse rates

% Find base temperatures and pressures
T_layers = To + [0, cumsum(L_layers.*diff(h_layers))];

P_layers = zeros(size(h_layers)); P_layers(1) = Po;
for jj = 2:numel(P_layers)
	Hb = h_layers(jj-1);
	Tb = T_layers(jj-1);
	Pb = P_layers(jj-1);
	L = L_layers(jj-1);
	H  = h_layers(jj);
	if L == 0
		P_layers(jj) = Pb*(exp(-GM_R*(H - Hb)/Tb));
	else
		P_layers(jj) = Pb*(T_layers(jj)/Tb)^(-GM_R/L);
	end
end

nz = numel(z);
T = zeros(nz,1);
P = zeros(nz,1);

for i = 1:nz

	if z(i) < 0
		warning('atmos:below_ground', ...
			'Input altitude must be greater than 0 (sea-level)');
		z(i) = 0.1;

	elseif z(i) >= 86000
		warning('atmos:above_atmos', 'Input altitude must be less than 84852m')
		z(i) = 86000;
	end


	h = E*z(i)/(E + z(i));	% Geopotential altitude

	% Find current layer
	dex = h >= h_layers;

	Tb = T_layers(dex);
	Pb = P_layers(dex);
	Hb = h_layers(dex);
	L  = L_layers(dex);
	
	T(i) = Tb + L*(h - Hb);

	if L == 0
		P(i) = Pb*(exp(-GM_R*(h - Hb)/Tb));
	else
		P(i) = Pb*(T(i)/Tb)^(-GM_R/L);
	end

end

T = reshape(T, size(z));
P = reshape(P, size(z));

% Density
if T <= 0
	error('Bah');
end
rho = P./(Ra*T);  % kg/m3

% Speed of Sound
a = sqrt(gamma*Ra*T);

% Viscosities
mu = mu_o*(T_o + C)./(T + C).*((T/T_o).^(3/2));
%mu = 1.458e-6*sqrt(T)/(1 + 110.4/T);

nu = mu./rho;

if show
	fprintf(1, '\nAt altitude %d m the atmospheric properties are:\n', z)
	fprintf(1, 'Temperature = %0.6g K\n', T)
	fprintf(1, 'Pressure = %0.6g kPa\n', P/1000)
	fprintf(1, 'Density = %0.6g kg/m3\n', rho)
	fprintf(1, 'Speed of sound = %0.6g m/s\n', a)
	fprintf(1, 'Dynamic viscosity = %0.6g Pa.s\n', mu)
	fprintf(1, 'Kinematic viscosity = %0.6g m^2/s\n\n', nu)
end

results = {T, P, rho, a, mu, nu};

nout = max(nargout,1)-1;

for i = 1:nout
    varargout(i) = results(i+1);
end

