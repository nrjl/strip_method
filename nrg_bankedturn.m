%--------------------------------------------------------------------------
%
% FUNCTION:		banked turn energy
%
% PURPOSE:		Determine the energy loss based on a trim method
%
% SYNTAX:		
%
% INPUTS:		V, phi, h
%
% OUTPUTS:		
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		trim_bankedturn
%--------------------------------------------------------------------------

global g
global PLANE_PARAM

g = 9.81;
alt = 100;
x_cg = -290*0.5;
y_cg = 0;
z_cg = 0;

% [SBXC_aero, SBXC_param] = SBXC_def(x_cg, y_cg, z_cg);

V_turn = 9:0.25:30;
nV = length(V_turn);

bank_range = 5:1:80;
nbank = length(bank_range);

edot = zeros(1, nV);
g_ratio = edot;
y_rate = edot;
m = PLANE_PARAM.m;

Xt = zeros(12, nV, nbank);
Ut = zeros( 4, nV, nbank);

failed = zeros(nbank, nV);

nrg_rate = zeros(nbank, nV);
glide_ratio = nrg_rate;
yaw_rate	= nrg_rate;

h_waitbar = waitbar(0, 'Calculating trimmed turn conditions');
for j = 1:nbank
	bank = bank_range(j);
	for i = 1:nV
		if i > 1
			if atan2(Xt(3,i-1,j), Xt(1,i-1,j)) < 5*pi/180
				[Xt(:,i,j), Ut(:,i,j), failed(j,i), Xdot]= trim_bankedturn(V_turn(i), bank,...
					alt, Xt(:,i-1,j), Ut(:,i-1,j));
			else
				[Xt(:,i,j), Ut(:,i,j), failed(j,i), Xdot]= trim_bankedturn(V_turn(i), bank,...
					 alt);
			end
		else
			[Xt(:,i,j), Ut(:,i,j), failed(j,i), Xdot] = trim_bankedturn(V_turn(i), bank, ...
				alt);
		end
		
		if any(isnan(Ut(:,i,j))) || failed(j,i)
			Xdot = [0 0 0 0 0 0 0 0 0 0 0 1]; failed(j,i) = 1; disp('Failed');
		end
		y_rate(i) = Xdot(9);
		edot(i) = Xdot(12)*m*g;
		g_ratio(i) = sqrt(sum(Xdot(10:11).^2))/Xdot(12);
		waitbar(((j-1)*nV + i)/(nV*nbank), h_waitbar);
	end
	
	nrg_rate(j,:)	= edot;
	yaw_rate(j,:)	= y_rate;
	glide_ratio(j,:)= g_ratio;
end

close(h_waitbar);
elevator	= squeeze(Ut(1,:,:));
aileron		= squeeze(Ut(2,:,:));
rudder		= squeeze(Ut(3,:,:));

%%

col = 'ybgrkmc';
sty = {':', '-', '--', '-.'};
m = PLANE_PARAM.m; g = 9.81;

colsty = cell(1, nbank);
leg_ent = cell(1, nbank);

figure(10); clf;
subplot(3,1,1); hold on; grid on
for i = 1:nbank
	colsty{i} = [col(mod(i, length(col))+1), sty{mod(i, length(sty))+1}];
	plot(V_turn, elevator(:,i), colsty{i});
	leg_ent{i} = sprintf('Bank angle = %0.3g deg', bank_range(i));
end
xlabel('Cruise velocity (m/s)'); ylabel('Elevator setting');	
legend(leg_ent); grid on;

subplot(3,1,2); hold on; grid on
for i = 1:nbank
	plot(V_turn, aileron(:,i), colsty{i});
end
xlabel('Cruise velocity (m/s)'); ylabel('Aileron setting'); grid on;

subplot(3,1,3); hold on; grid on
for i = 1:nbank
	plot(V_turn, rudder(:,i), colsty{i});
end
xlabel('Cruise velocity (m/s)'); ylabel('Rudder setting');  grid on;

figure(11); clf;
subplot(2,1,1); hold on; grid on
for i = 1:nbank
	plot(V_turn, nrg_rate(i,:), colsty{i});
end
xlabel('Cruise velocity (m/s)'); ylabel('Energy loss rate (W)');
legend(leg_ent); grid on;

subplot(2,1,2); hold on; grid on
for i = 1:nbank
	plot(V_turn, glide_ratio(i,:), colsty{i});
end
xlabel('Cruise velocity (m/s)'); ylabel('Glide ratio'); grid on;

%%
% Yaw rate versus bank angle
figure(12); clf; hold on; grid on

if nV <= 10
	
	vlegend = cell(1, nV);
	for i = 1:nV
		lsty = [col(mod(i, length(col))+1), sty{mod(i, length(sty))+1}];
		plot(bank_range, yaw_rate(:,i)'*180/pi, lsty);
		vlegend{i} = sprintf('V = %0.4g m/s', V_turn(i));
	end
	xlabel('Bank angle (deg)'); ylabel('Yaw rate (deg/s)');
	legend(vlegend); grid on;
else
	h_yaw = surf(V_turn, bank_range, yaw_rate*180/pi, 'MeshStyle', 'column');
	xlabel('Airspeed(m/s)'); ylabel('Bank angle (deg)');
	zlabel('Yaw rate (deg/s)');
end
	

% Turn efficiency
turn_eff = yaw_rate*m*g*180/pi./nrg_rate;

figure(14); clf;
h_eff = surf(V_turn, bank_range, turn_eff, 'MeshStyle', 'column'); %shading('interp')
colormap hot;
xlabel('Airspeed (m/s)'); ylabel('Bank angle (deg)');
zlabel('Turn efficiency (deg yaw/m altitude lost)');