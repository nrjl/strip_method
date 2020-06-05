%% Initialisation & globals
initialise_SBXC

global g psi_w
g = 9.81;
psi_w = rectify(pi);
wind_dir = psi_w;

%% Domain
h_start = 100;
V_start = 10:1:20;
nV = length(V_start);
gamma_a = (-50:2.5:50)*pi/180;
ng = length(gamma_a);

dt = 0.1;

%% Variable allocation
X0 = zeros(12, nV);
U0 = zeros(4, nV);
V_end = zeros(ng, nV);
h_end = zeros(ng, nV);

for j = 1:nV
	[X0(:,j), U0(:,j)] = trim_controls(V_start(j), h_start);
	
	for i = 1:ng
		X0(8, j) = gamma_a(i) + atan2(X0(3,j), X0(1,j));

		dt_block = dt/5;
		X_block = euler_stupid(X0(:,j), U0(:,j), dt, dt_block);

		V_end(i,j) = sqrt(sum(X_block(1:3,end).^2, 1));
		h_end(i,j) = -X_block(12,end);
		fprintf('.');
	end
	fprintf('\n');
end

m = PLANE_PARAM.m;
Power = (.5*m*V_end.^2 + m*g*h_end - .5*m*(ones(ng,1)*V_start).^2 - m*g*h_start)/dt;
Accel = (V_end - (ones(ng,1)*V_start))./dt;

%% Plotting
figure(1); clf;

subplot(1,2,1);
surf(V_start, gamma_a*180/pi, Power);
xlabel('Airspeed, {\itV} (m/s)'); ylabel('Climb Angle, {\it\gammai} (deg)'); ...
	zlabel('Power, {\itdE/dt} (J/s, W)');
title('Power in Steady Climb');

subplot(1,2,2);
surf(V_start, gamma_a*180/pi, Accel);
xlabel('Airspeed, {\itV} (m/s)'); ylabel('Climb Angle, {\it\gammai} (deg)'); ...
	zlabel('Air acceleration, {\itdV/dt} (m/s/s)');
title('Air acceleration in Steady Climb');

%% Overplotting
% colour = @(dex) [0, 255, 0]./255 + [1,-1,0]*dex +
% [0,0,-2*abs(dex-0.5)+1];

VL = find(V(1,:) >= V_start(1)); VL = VL(1);
VH = find(V(1,:) <= V_start(end)); VH = VH(end);

gL = find(gammai(:,1) >= gamma_a(1)); gL = gL(1);
gH = find(gammai(:,1) <= gamma_a(end)); gH = gH(end);

figure(1); subplot(1,2,1); hold on;
h2 = surf(V(1,VL:VH), gammai(gL:gH, 1)*180/pi, dE_dt(gL:gH, VL:VH));
set(h2, 'facealpha', 0.5, 'edgealpha', 0.5);

subplot(1,2,2); hold on;
h3 = surf(V(1,VL:VH), gammai(gL:gH, 1)*180/pi, dV_dt(gL:gH, VL:VH));
set(h2, 'facealpha', 0.5, 'edgealpha', 0.5);