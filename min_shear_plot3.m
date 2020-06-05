%% Script: minimum shear 3d plot

%% SETUP
% if ~exist('SBXC_aero', 'var')
% 	x_cg = -0.5*290;
% 	y_cg = 0;
% 	z_cg = 0;
% 	[SBXC_aero, SBXC_param] = SBXC_def(x_cg, y_cg, z_cg);
% end

%% MAIN
alt = 100;
vel_range = 10:0.05:20;				% Airspeed range
gamma_range = (5:0.5:30)*pi/180;	% Climb angle range

nV = length(vel_range);
ng = length(gamma_range);

D3 = zeros(1, nV);

beta3 = zeros(nV, ng);
Edot_wind3 = zeros(nV, ng);

gamma = zeros(1, nV);
beta = zeros(1, nV);
Edot_wind = zeros(1, nV);
D = zeros(1, nV);

for i = 1:nV	
	[gamma_range, beta3(i,:), Edot_wind3(i,:), D3(i)] = ...
		min_shear(vel_range(i), alt, gamma_range);	
	[gamma(i), beta(i), Edot_wind(i), D(i)] = ...
		min_shear(vel_range(i), alt);
end

%% Plotting
figure(13); clf;

C_map = beta3' - (beta3(:,1)*ones(1, size(beta3, 2)))';
beta_surf = surf(vel_range, gamma_range*180/pi, beta3', C_map);
set(beta_surf, 'EdgeColor', 'none')
xlabel('Airspeed (m/s)'); ylabel('Climb angle, \gamma (deg)'); 
zlabel('Required wind shear, \beta (1/s)');

hold on; plot3(vel_range, gamma*180/pi, beta, '-k')

% figure(14); clf;
% 
% E_surf = surf(vel_range, gamma_range*180/pi, Edot_wind3');
% set(E_surf, 'EdgeColor', 'none')
% xlabel('Airspeed (m/s)'); ylabel('Climb angle, \gamma (deg)'); 
% zlabel('Wind power (W)');
