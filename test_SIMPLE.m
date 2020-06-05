%% Test Scratchpad
addpath graphics wind_profiles wing_data SBXC static_functions
initialise_SBXC

% reset(RandStream.getDefaultStream);
% defaultStream = RandStream.getDefaultStream;
defaultStream = RandStream.getGlobalStream;
savedState = defaultStream.State;

global g psi_w
g = 9.81;
psi_w = rectify(pi);
wind_dir = psi_w;

oplot = 0;	% overplot flag
colrs = 'bgrmy'; lss = {'-', '--', '-.', ':', '-'};
col = colrs(oplot+1);
ls = lss{oplot+1};

start_time = cputime;

alt = 100;
V_initial = 20;
bank_initial = 0;

V_target = V_initial;
bank_target = bank_initial;
V_cycle = V_target;
KE_cycle = 0.5*PLANE_PARAM.m*V_cycle^2;

if bank_target ~= 0
	[X0, U0] = trim_bankedturn(V_target, bank_target, alt);
	bank_target = bank_target*pi/180;
else
	[X0, U0] = trim_controls(V_target, alt);
end

% Flightmodes:
%	0 straight and level
flightmode = 0; alt_target = alt;
travel_dir = 0;
travel_vector = [cos(travel_dir); sin(travel_dir)];
lookahead_distance = 500;

going_down = false;
max_alt = alt; min_alt = alt;
gam = 0*pi/180;
high_travel = false;		% Whether the travel is done during the high turn (true) or low turn (false)

banked_turn = 55*pi/180;
turn_cost = 150; % (J/rad)

flags = [1, -1];
section_number = 2;
% X0(8) = gam + atan2(X0(3), X0(1));

%% Load beta shear variables
if ~exist('vel_range', 'var')
	load energy_gain.mat
end

%%
figure(1);
if oplot
	hold on
	delete(h_sbxc);
else
	clf
	[XX, YY, ZZ] = meshgrid(linspace(60, 60, 1), ...
							linspace(-50, 50, 10), ...
							linspace(-130, -90, 10));
	V_w = wind_field([XX(:)'; YY(:)'; ZZ(:)']);
	U_wind = reshape(V_w(1,:), size(XX));
	V_wind = reshape(V_w(2,:), size(YY));
	W_wind = reshape(V_w(3,:), size(ZZ));
	vscale = PLANE_PARAM.fuse_l/V_target*2;
	h_quiver = quiver3(XX, -YY, -ZZ, U_wind*vscale, V_wind*vscale, W_wind*vscale, 0);
	set(h_quiver, 'color', [0.6 0 0]);
end

%%
set(1, 'Name', 'Trajectory')
XX0 = dim_fix(X0);
h_sbxc = SBXC_handle(XX0, SBXC_aero, 1);
set(gca, 'Xdir', 'normal', 'YDir', 'normal', 'ZDir', 'normal')
ax_all =[XX0(10)-10, XX0(10)+10, XX0(11)-10, XX0(11)+10, XX0(12)-10, XX0(12)+10];
axis(ax_all)


%% Control stuff
servo_rate = 4;  % 2
turn_dir = 1;

% Velocity Control
Kp_V = 	0.2;
Ki_V = 	0.07;
Kd_V = 	0.2;
V_int = 0;

% Climb Control
Kpp = 	-4;
Kip = 	-4;
Kdp = 	0.2;
climb_int = 0;

% Bank angle control
Kp_bank =	8;
Ki_bank =	1.0;
Kd_bank =	1.0;
bank_int = 0;

% Sideslip control
Kp_slip = -0.8;
Ki_slip = -0.25;
Kd_slip = -0.03;
slip_target = 0;
slip_int = 0;

% Heading control
Kp_hb = 1.2;
Ki_hb = 0;
Kd_hb = 0.1;
Kp_hs = 1.0;
Ki_hs = 0;
Kd_hs = 0.25;
head_int = 0;

% Altitude control
Kp_alt = 0.8*pi/180;

% Bank angle parameters
K_bank = [80.42674848184565, 1.2178304290749726];
V_min = 10;

% Sensor errors
sigma_p = 0.3;			% Pitot static velocity error std (m/s)
sigma_a = 0.2*pi/180;	% Alpha vane error std (rad)
sigma_b = 0.2*pi/180;	% Beta vane error std (rad)
sigma_gps = 0.0;		% GPS velocity error

%% Integration Loop
dt = 0.1;
t0 = 0;
tf = 20;

np = round((tf - t0)/(dt/2^4) + 1);
nc = round((tf - t0)/dt + 1);

% Pre-allocation
U_full = U0*ones(1, nc);

X_mov  =	zeros(length(X0), nc);
X_full =	zeros(length(X0), np);
uvw_wind =	zeros(3, np);

t_full =	zeros(1, np);
V =			zeros(1, np);

% Control variables initialisation
q_err =		zeros(1, nc);
V_err =		zeros(1, nc);
climb_err = zeros(1, nc);
bank_err =	zeros(1, nc);
slip_err =	zeros(1, nc);
roll_err =	zeros(1, nc);
head_err =	zeros(1, nc);
slip =		zeros(1, nc);

modetimes = zeros(1, 5);
moderow = 1;

% Sensor variables init
pitot =		zeros(1, nc);
alpha_vane=	zeros(1, nc);
beta_vane =	zeros(1, nc);
obs_z = [];
obs_wind = [];
n_samples = 100;

Cr = eye(4); Ct = eye(4);

bar_handle = waitbar(0, 'Integrating flight, please wait...');

% Initial conditions
ctime = t0; modetimes(1,1) = t0;
X_full(:,1) = X0;
X_mov(:,1)  = X_full(:,1);
uvw_wind(:,1) = X0(1:3) - calc_Ceb(X0(7:9))'*wind_field(X0(10:12));

t_full(1) = t0;
V(1) = sqrt(sum(uvw_wind(1:3,1).^2));
slip(1) = asin(uvw_wind(2,1:1)/V(1));

wind_est = [1; pi/2; -pi+5*pi/180];
current_est = wind_est;

plot_skip = fix(2.5/dt);
wskip = ceil((tf - t0)/dt/100);

istart = 2;
lvl = 1;
rk4tol = 1e-3;

% Control initialisation
V_on	= true;	V_target = V_initial;
climb_on= false;
climba_on=false;		climb_target = 0;%gam;
alt_on	= false;
bank_on	= true;
slip_on	= false;
head_on = true;		head_target = pi/2;

for i = 2:(tf-t0)/dt+1

% 	if i == 2 || lvl <= 1
% 		[X_block, lvl] = rk4_variable(X_full(:,istart-1), U_full(:,i-1), ...
% 			dt, rk4tol, [], flags(1), flags(2)); %0,0 1, -1
% 		dt_block = dt/(2^lvl);
% 	else
% 		ctol = rk4tol*(0.75^(lvl-1));
% 		[X_block, lvl] = rk4_variable(X_full(:,istart-1), U_full(:,i-1), ...
% 			dt/(2^(lvl-1)), ctol, [], lvl, 1);
% 		dt_block = dt/(2^lvl);
% 	end

	dt_block = dt/5;
	X_block = euler_stupid(X_full(:,istart-1), 	U_full(:,i-1), dt, dt_block);
	
	
	istop = istart + dt/dt_block - 1;
	
	X_full(:,istart:istop) = X_block;
	X_mov(:,i) = X_block(:,end);
	X_full(7:9,istart:istop) = rectify(X_full(7:9,istart:istop));
	t_full(istart:istop) = (dt_block:dt_block:dt) + t_full(istart-1);
	Ceb = calc_Ceb(X_full(7:9,istop));
		
	for k = istart:istop
		uvw_wind(:,k) = X_full(1:3,k) - ...
			calc_Ceb(X_full(7:9,k))'*wind_field(X_full(10:12,k));
	end	
	
	V(istart:istop) = sqrt(sum(uvw_wind(:,istart:istop).^2, 1));
	v_earth = Ceb*X_full(1:3,istop);
	slip(i) = asin(uvw_wind(2,istop)./V(istop));
	
	track_bearing = atan2(X_full(11,istop)-X0(11), X_full(10, istop) - X0(10));
	intrack_distance = dot(X_full(10:11,istop)-X0(10:11), travel_vector);
	target_point = travel_vector*(intrack_distance + lookahead_distance);
	target_bearing = atan2(target_point(2)-X_full(11,istop), ...
		target_point(1)-X_full(10, istop));

	
	% ---------- CONTROL STUFF ---------- %	
	
    g_temp = gam;
    
    V_on	= false;
    climb_on= true;		climb_target = g_temp;
    alt_on	= false;
    bank_on	= true;		bank_target = 0;
    slip_on	= false;
    head_on = true;		head_target = rectify(wind_dir - pi);%	0; %

	
	%% STALL PROTECTION
	if V(istop) < 9
		fprintf(1, 'Stall protection activated\n');
		V_on = true; V_target = 12;
		climb_on = false;
	end
	
	% Airspeed control (elevator actuated)
	if V_on
		V_err(i) = V_target - V(istop);
		
		if abs(V_err(i)) < 5
			climba_on = 0;
			V_int = V_int + V_err(i)*dt*(abs(U_full(1,i-1))~=1);

			U_full(1,i) = U_full(1,1) + ...
				Kp_V*V_err(i) + ...
				Ki_V*V_int + ...
				Kd_V*(V_err(i) - V_err(i-1))/dt;

		else
			climba_on = 1;
			climb_target = -sign(V_err(i))*30*pi/180;
		end
	end
	
	% Altitude control (climb angle controlled)
	if alt_on; climb_target = (alt_target+X_full(12,istop))*Kp_alt; end
	
	% Climb control (elevator actuated) (climba is air-relative, climb is
	% inertial climb wrt airspeed)
	if climb_on || climba_on
		if climb_on
			c_climb = asin(-v_earth(3)/V(istop));
		else
			c_climb = X_full(8,istop) - atan2(uvw_wind(3,istop), uvw_wind(1,istop));
		end
		climb_err(i) = climb_target - c_climb;
		climb_int = climb_int + climb_err(i)*dt*(abs(U_full(1,i-1))~=1);
		U_full(1,i) =	U_full(1,1) + ...
						Kpp*climb_err(i) + ...
						Kip*climb_int + ...
						Kdp*X_full(5,istop);
	end

	
	% Heading control (bank angle and rudder actuated)
	if head_on
		slip_on = false; bank_on = true;
		head_err(i) = rectify(head_target - X_full(9,istop));
		head_int = head_int + head_err(i)*dt;
		
		if abs(head_err(i)) > pi/3
			bank_target = sign(head_err(i))*(bank_exp(K_bank, V_min, V(istop))-(V(istop)-V_cycle)*pi/180);
			slip_on = true; slip_target = 0;
		else			
			bank_target = 	Kp_hb*head_err(i) + ...
				Ki_hb*head_int + ...
				Kd_hb*(head_err(i) - head_err(i-1))/dt;
			bank_target = saturate(bank_target, banked_turn);
			U_full(3,i) =	U_full(3,1) +...
				Kp_hs*head_err(i) + ...
				Ki_hs*head_int + ...
				Kd_hs*(head_err(i) - head_err(i-1))/dt;
		end
	end
	
	% Bank angle control (aileron actuated)	
	if bank_on && i >= 3
		bank = X_full(7, istop);
		bank_err(i) = rectify(bank_target - bank);
		bank_int = bank_int + bank_err(i)*dt*(abs(U_full(2,i-1))~=1);
		U_full(2,i) =	U_full(2,i-1) + ...
						Kp_bank*(bank_err(i)-bank_err(i-1)) + ...
						Ki_bank*bank_err(i) + ...
						Kd_bank*(bank_err(i) - 2*bank_err(i-1) + bank_err(i-2));
	end
	
	% Sideslip control (rudder actuated)
	if slip_on
		slip_err(i) = rectify(slip_target - slip(i));
		slip_int = slip_int + slip_err(i)*dt*(abs(U_full(3,i-1))~=1);
		U_full(3,i) = U_full(3,1) + (Kp_slip*slip_err(i) + Ki_slip*slip_int + ...
			Kd_slip*(slip_err(i) - slip_err(i-1))/dt);
	end
	
	% Saturate limits
	U_full(:,i) = saturate(U_full(:,i), 1);	
	U_full(:,i) = servo_check(U_full(:,[i-1,i]), servo_rate, dt);
	
	
	% --------------- PLOTTING --------------- %	
	% Plot path
	set(0,'CurrentFigure',1);
	XXnow = dim_fix(X_full(:,istart-1:istop));
	plot3(XXnow(10, :), XXnow(11, :), XXnow(12, :), [col, '-']);
	ax_low = min([ax_all([1,3,5])', XXnow(10:12,end)], [], 2);
	ax_high= max([ax_all([2,4,6])', XXnow(10:12,end)], [], 2);
	ax_all = zeros(1,6);
	ax_all([1,3,5]) = ax_low'; ax_all([2,4,6]) = ax_high';
	axis(ax_all);

	Cr = eye(4);
	Cr(1:3, 1:3) = calc_Ceb(XXnow(7:9,end));% Rotation transformation
	Ct = eye(4);
	Ct(1:3,4) = XXnow(10:12, end);			% Translate to actual position
	set(h_sbxc, 'Matrix', Ct*Cr);
	
	if ~mod(i, wskip)
		wait_str = sprintf('%1.0f%% complete', i/nc*100);
		waitbar(i/nc, bar_handle, wait_str);
	end
	
	if ~mod(i, 1/dt)
		plot3(XXnow(10,end), XXnow(12,end), XXnow(12,end), '.')
		if ~mod(i, 5/dt)
			text(XXnow(10,end), XXnow(11,end), XXnow(12,end),...
				sprintf('%0.2gs', i*dt+t0), 'FontSize', 6);
		end
	end
	
	if moderow > 5 && flightmode == 0; break; end
	
	drawnow;	
	istart  = istop+1;
end

%%
close(bar_handle);

total_time = cputime - start_time;
fprintf(1, 'Simulation complete in %0.5f sec (%0.4f�real time)\n\n', ...
	total_time, total_time/(t_full(istop) - t_full(1)));

%% Cut matrices for unfinished sims
X_full		= X_full(:,1:istop);
uvw_wind	= uvw_wind(:,1:istop);
t_full		= t_full(:,1:istop);
V			= V(:,1:istop);
tv			= t0:dt:round(t_full(end)/dt)*dt;
U_full		= U_full(:, 1:length(tv));
X_mov		= X_mov(:,1:length(tv));

%% Plotting
plot_states(X_full, U_full, t_full, tv, PLANE_PARAM.m, oplot, col, ls, modetimes)

%% Movie
figure(9); clf;
fskip = 2; 
XYZ = {XX, YY, ZZ};
UVW = {U_wind*vscale, V_wind*vscale, W_wind*vscale};

[M] = plane_movie(SBXC_aero, fskip, X_mov, U_full, col, XYZ, UVW, 0.6);

