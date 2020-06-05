%%
deriv_surf = diff(turn_eff, 1, 1);
deriv_surf = deriv_surf.*(abs(deriv_surf) < 2.2);

dderiv_surf = diff(deriv_surf, 1, 1);
dderiv_surf = dderiv_surf.*(abs(dderiv_surf) < 0.2);

m = SBXC_param.m;

%%
eff_lim = zeros(3, nV);
bank_max= zeros(1, nV);
eff_max = zeros(1, nV);

figure(1); clf; hold on;

for i = 1:nV
	
	eff_lim(1,i) = turn_eff(1,i);
	
	for j = 3:nbank
		if turn_eff(j,i) <= 0
			eff_lim(2,i) = turn_eff(j-1,i)*((nbank-j-2)/(nbank-1)*0.5+0.5);
			eff_lim(3,i) = j-1;
			break
		elseif j == nbank
			eff_lim(2,i) = turn_eff(j,i)*((nbank-j-2)/(nbank-1)*0.5+0.5);
			eff_lim(3,i) = j;
		end
	end
% 	d = turn_eff(1:eff_lim(3,i),i)' - linspace(eff_lim(1,i), eff_lim(2,i), eff_lim(3,i));
	v_hat = 1/sqrt((eff_lim(2,i) - eff_lim(1,i))^2 + (bank_range(eff_lim(3,i)) - bank_range(1))^2);
	v_hat = v_hat * [-(eff_lim(2,i) - eff_lim(1,i)); (bank_range(eff_lim(3,i))-bank_range(1))];
	d = [bank_range(1:eff_lim(3,i))', turn_eff(1:eff_lim(3,i), i)]*v_hat;
	max_index = find( d == max(d));
	if any(d) ~= 0; max_index = max_index(end);	else max_index = max_index(1); end
	bank_max(i)  = bank_range(max_index(1));
	eff_max(i)   = turn_eff(max_index(1), i);
	plot(bank_range, turn_eff(:, i), 'b-', [bank_range(1), bank_range(eff_lim(3,i))], [eff_lim(1,i), eff_lim(2,i)], 'r-',  bank_range(max_index(1)), eff_max(i), 'ko')
end

%%
figure(14); hold on;
% if exist('h_maxeff', 'var')
% 	delete(h_maxeff)
% 	clear h_maxeff
% end
h_maxeff = plot3(V_turn, bank_max, eff_max, 'Color', [0 0.8 0]);
xlabel('Airspeed (m/s)'); ylabel('Bank angle (deg)'); zlabel('Turn efficiency (deg/m)');

%%
figure(15); clf;
h_surf2 = surf(V_turn, bank_range, turn_eff, 'MeshStyle', 'column');
set(h_surf2, 'CData', [deriv_surf(1,:); deriv_surf]);
colormap jet
hold on; h_maxeff2 = plot3(V_turn, bank_max, eff_max, 'Color', [0 0 0.8]);
xlabel('Airspeed (m/s)'); ylabel('Bank angle (deg)'); zlabel('Turn efficiency (deg/m)');

%%
figure(16); clf;
h_surf3 = surf(V_turn, bank_range, turn_eff, 'MeshStyle', 'column');
set(h_surf3, 'CData', [dderiv_surf(1,:); dderiv_surf; dderiv_surf(end, :)]);
colormap jet
hold on; h_maxeff3 = plot3(V_turn, bank_max, eff_max, 'Color', [0 0 0.8]);
xlabel('Airspeed (m/s)'); ylabel('Bank angle (deg)'); zlabel('Turn efficiency (deg/m)');

%%
non_min = find(bank_max ~= min(bank_range));
figure(17); clf; plot(V_turn(non_min), bank_max(non_min), 'ro'); hold on;

if ~exist('V_test', 'var')
	V_test = linspace(min(V_turn(non_min)), max(V_turn(non_min)), 101);
	
	% P4_2 = polyfit(V_turn(18:end), bank_max(18:end), 4);
	% plot(V_test, polyval(P4_2, V_test), 'Color', [0 0.8 0], 'LineStyle',
	% '--');
	
	K = [80.42674848184565, 1.2178304290749726];
	bank_out = bank_exp(K, 10, V_test)*180/pi;
	
	turn_eff_fit = zeros(1, length(V_test));

	for i = 1:length(V_test)
		[X, U, ff, Xdot] = trim_bankedturn(V_test(i), bank_out(i), 100);
		turn_eff_fit(i) = Xdot(9)/Xdot(12)*180/pi;
	end
end

plot(V_test, bank_out, 'b-');
xlabel('Airspeed (m/s)'); ylabel('Bank angle (deg)');

figure(14); hold on;
plot3(V_test, bank_out, turn_eff_fit, 'k-', 'LineWidth', 2.0);

figure(15); hold on;
plot3(V_test, bank_out, turn_eff_fit, 'k-', 'LineWidth', 2.0);