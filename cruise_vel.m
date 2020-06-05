%--------------------------------------------------------------------------
%
% FUNCTION:		cruise_vel
%
% PURPOSE:		Determine the energy loss based on a trim method
%
% SYNTAX:		
%
% INPUTS:		h	- Altitude (m)
%
% OUTPUTS:		
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		trim_controls
%--------------------------------------------------------------------------
global g
global PLANE_AERO PLANE_PARAM
g = 9.81;
h = 100;
[PLANE_AERO, PLANE_PARAM] = SBXC_def(-0.5*290, 0, 0);

V_cruise = [10.5:0.1:11.3, 11.35:0.05:12.3, 12.4:0.1:14, 14.5:0.5:20]; %[9.9:0.1:11.5, 11.55:0.05:12.3, 12.4:0.1:14];
np = length(V_cruise);
edot = zeros(1, np);
g_ratio = edot;

cg_range = 0.5; %[0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7];
mass_range = (12)*0.45359237;
ncg = length(mass_range);

Xt = zeros(12, np, ncg);
Ut = zeros( 4, np, ncg);

y_cg = 0;
z_cg = 0;

nrg_rate = zeros(ncg, np);
elevator = nrg_rate;
glide_ratio = nrg_rate;
leg_ent = cell(1, ncg);

h_waitbar = waitbar(0, 'Calculating trimmed cruise condition');

for j = 1:ncg
% 	x_cg = -cg_range(j)*290;
% 	[PLANE_AERO, PLANE_PARAM] = SBXC_def(x_cg, y_cg, z_cg);
	PLANE_PARAM.m = mass_range(j);

	for i = 1:np
		[Xt(:,i,j), Ut(:,i,j)] = trim_controls(V_cruise(i), h);
		Xdot = state_rates(Xt(:,i,j), Ut(:,i,j));
		edot(i) = Xdot(12)*PLANE_PARAM.m*g;
		g_ratio(i) = Xdot(10)/Xdot(12);
	end
	
	nrg_rate(j,:) = edot;
	elevator(j,:) = Ut(1,:,j);
	glide_ratio(j,:) = g_ratio;
	leg_ent{j} = sprintf('mass = %0.3gkg', mass_range(j));
	waitbar(j/ncg, h_waitbar);
end
close(h_waitbar);

%% Theoretical results
e = 0.85;
[S, AR] = plane_properties(PLANE_AERO, PLANE_PARAM);
[T, P, rho] = atmos(h);
% Cd0_p = polyfit([12 15 20], [0.011902 .011214 .010247], 0);
% polyval(Cd0_p, V_cruise)

pre = pi*AR*e;
q = 0.5*rho*V_cruise.^2;
gamma = asin(1./(2*mass_range*g)*(-pre*q*S + sqrt(...
	pre^2*q.^2*S^2 + 4*pre*q.^2*S^2.*SBXC_param.Cd0 + 4*mass_range.^2*g^2)));


%% Plotting
col = 'ybgrkmc';
sty = {':', '-', '-.', '--'};
mark = 'ox+*sdv^<>ph';

figure(10); clf; hold on
colsty = cell(1, ncg);
for i = 1:ncg
	colsty{i} = [col(mod(i, length(col))+1), ...
		mark(mod(i, length(mark))+1), sty{mod(i, length(sty))+1}];
	plot(V_cruise, elevator(i,:), colsty{i}, 'MarkerSize', 3);
end
xlabel('Cruise velocity (m/s)'); ylabel('Elevator setting');
title('Elevator position for trimmed cruise condition')
legend(leg_ent); grid on;

%%
figure(11); clf;
subplot(2,1,1); hold on
for i = 1:ncg
	plot(V_cruise, nrg_rate(i,:)./mass_range(i)/g, colsty{i}, 'MarkerSize', 3);
end
plot(V_cruise, V_cruise.*sin(gamma), 'r--');
leg_ent{j+1}  = 'Theoretical Prediction';
xlabel('Cruise velocity (m/s)'); ylabel('Sink rate (m/s)');
title('Sink rate polar');
legend(leg_ent); grid on;

% figure(12); clf; hold on
subplot(2,1,2); hold on
for i = 1:ncg
	plot(V_cruise, glide_ratio(i,:), colsty{i}, 'MarkerSize', 3);
end
plot(V_cruise, 1./tan(gamma), 'r--');
leg_ent{j+1}  = 'Theoretical Prediction';
xlabel('Cruise velocity (m/s)');
ylabel('Glide ratio');
legend(leg_ent); grid on;