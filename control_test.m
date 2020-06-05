% clear all
global g
global PLANE_AERO PLANE_PARAM
g = 9.81;

clear bank_target

run_sim = 0;
oplot = 0;
col = 'g'; ls = '--';

x_cg = -0.5*290;
y_cg = 0;
z_cg = 0;

[SBXC_aero, SBXC_param] = SBXC_def(x_cg, y_cg, z_cg);
plane_properties(SBXC_aero, SBXC_param)
PLANE_AERO = SBXC_aero;
PLANE_PARAM = SBXC_param;

alt = 100;
V_target = 12;
% bank_target = 0;

if exist('bank_target', 'var')
	[X0, U0] = trim_bankedturn(V_target, bank_target, alt);
else
	bank_target = 0;
	[X0, U0] = trim_controls(V_target, alt);
end

[A, B] = state_space(@(X, U) state_rates(X, U), X0, U0);

% u w q theta x z
Alon_rows = [1,3,5,8,10,12];
Blon_rows = [1,4];
Alon = A(Alon_rows, Alon_rows);
Blon = B(Alon_rows, Blon_rows);

% v p r phi psi y
Alat_rows = [2,4,6,7,9,11];
Blat_rows = [2,3];
Alat = A(Alat_rows, Alat_rows);
Blat = B(Alat_rows, Blat_rows);

disp('LATERAL')
damp(Alat)
disp('LONGITUDINAL')
damp(Alon)

lat_eig = eig(Alat);
lon_eig = eig(Alon);

%% ----- TRANSFER FUNCTIONS ----- %%

% Longitudinal
C = zeros(1, size(Alon,1));
D = zeros(1, size(Blon,2));

% Elevator effect on airspeed
iu = 1;					% Elevator, de
C(1) = 1;               % Forward speed, u

[num_de_u, den_de_u] = ss2tf(Alon, Blon, C, D, iu);
[Z_de_u, P_de_u, K_de_u] = ss2zp(Alon, Blon, C, D, iu);

tf_de_u = zpk(Z_de_u, P_de_u, K_de_u);

% Lateral
C = zeros(1, size(Alat,1));
D = zeros(1, size(Blat,2));

% Aileron effect on roll rate
iu = 1;					% Aileron, da
C(2) = 1;               % Roll rate, p

[num_da_p, den_da_p] = ss2tf(Alat, Blat, C, D, iu);
[Z_da_p, P_da_p, K_da_p] = ss2zp(Alat, Blat, C, D, iu);

tf_da_p = zpk(Z_da_p, P_da_p, K_da_p)


%% ----- INTEGRATE SS APPROXIMATION ----- %%
if run_sim == 1
	dt = 0.1;
	t0 = 0;
	tf = 20;
	np = (tf-t0)/dt + 1;

	X = zeros(length(X0), np); X(:,1) = X0;
	uvw_wind = zeros(3, np);
	U = U0*ones(1, np);
	
	V = zeros(1,np);		V_int = 0;		V_err = zeros(1,np);
	slip = zeros(1,np);		slip_int = 0;
	bank = zeros(1,np);		bank_int = 0;

	Xdot0 = state_rates(X0, U0);

	V_target = 13;
	bank_target = 0;
	slip_target = 0;
% 	U(1, 5:8) = 0.1*ones(1,4);

	for i = 2:np
		Xdot = A*(X(:,i-1)-X0) + B*(U(:,i-1)-U0);
		Xdot(10:12) = Xdot(10:12) + Xdot0(10:12);
		X(:,i) = X(:,i-1) + Xdot*dt;

		uvw_wind(:,i) = X(1:3,i) - ...
			calc_Ceb(X(7:9,i))'*wind_field(X(10:12,i));

	% Airspeed control (elevator actuated)
	if V_on
		V_err(i) = V_target - V(i);
		V_int = V_int + V_err(i)*dt*(abs(U(1,i-1))~=1);
		
		U(1,i) = U(1,1) + ...
			Kp_V*V_err(i) + ...
			Ki_V*V_int + ...
			Kd_V*(V_err(i) - V_err(i-1))/dt;
	end

		U(:,i) = servo_check(U(:,[i-1,i]), servo_rate, dt);		
	end
	disp('Linear sim complete');

	%% Plotting

	t = t0:dt:tf;
	m = SBXC_param.m;

	plot_states(X, U, t, t, m, oplot, col, ls)
end