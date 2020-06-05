% Compare theory and sim for banked turns
plot_polars;

lift_and_drag;

load dat\turn_nrg_big_edit.mat
m = PLANE_PARAM.m;

%%
pre = pi*AR*e;
q = 0.5*rho*V_turn.^2;
phi = bank_range*pi/180;

Cd0_p = polyfit([12 15 20], [0.011902 .011214 .010247], 0);

gamma = asin(1/(2*m*g)*(-pre*cos(phi').^2*q*S + sqrt(...
	pre^2*cos(phi').^4*q.^2*S^2 + 4*pre*cos(phi').^2*(q.^2*S^2.*polyval(Cd0_p, V_turn)) + 4*m^2*g^2)));
CL_all = m*g*cos(gamma)./(cos(phi')*q*S);

% Theoretical glide ratio
tGR = 1./tan(gamma); %+ polyval(p_Cd, CL));

figure(8); clf; hold on;
h_tGR = surf(V_turn, bank_range, tGR, tGR-max(glide_ratio(:))/2, ...
	'FaceAlpha', 0.4, 'EdgeAlpha', 0.4, 'MeshStyle', 'row');
glide_ratio(find(glide_ratio == 0)) = NaN;
h_sGR = surf(V_turn, bank_range, glide_ratio, 'MeshStyle', 'row', 'EdgeAlpha', 0.4);
xlabel('Airspeed(m/s)'); ylabel('Bank angle (deg)');
zlabel('Glide ratio');

figure(9); clf; hold on;
yaw_theory = CL_all*0.5*rho*S.*(sin(phi')*V_turn)/m;
h_tYR = surf(V_turn, bank_range, yaw_theory*180/pi, 180/pi*...
	(yaw_theory./3), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
h_sYR = surf(V_turn, bank_range, yaw_rate*180/pi, 'MeshStyle', 'column', 'EdgeAlpha', 0.4);
xlabel('Airspeed(m/s)'); ylabel('Bank angle (deg)');
zlabel('Yaw rate (deg/s)');

figure(10); clf; hold on;
theory_eff = (tan(phi')*(1./(V_turn.^2))*g.*tGR*180/pi);
turn_eff(find(turn_eff == 0)) = NaN;
h_tTE = surf(V_turn, bank_range, theory_eff, theory_eff./2+max(turn_eff(:)), 'MeshStyle', 'column', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
h_sTE = surf(V_turn, bank_range, turn_eff, 'MeshStyle', 'column');
xlabel('Airspeed(m/s)'); ylabel('Bank angle (deg)');
zlabel('Turn efficiency (deg/m)');
