function [min_alt, max_alt, direction] = find_ds_limits(obs_loc, obs_wind, beta_required)
% Function to calculate limits for dynamic soaring

if isempty(obs_wind)
	direction = pi;
else
	direction = atan2(mean(sin(obs_wind(3,:))), mean(cos(obs_wind(3,:))));
end

%% GP input variables
X = obs_loc;			% training points
y = obs_wind(1,:)';     % training outputs
sigma_f = 1.0;			% process std
sigma_n = 0.5;			% observation std
l = 2;                  % length scale
n_obs = numel(y);       % Number of observations (length y)

%% Test Points
n_test = 100;            % Number of target points
x_star = linspace(min(X), max(X), n_test);

%%
%----- COVARIANCE FUNCTION -----%
K_maker = @(x1,x2, sigma_f)  sigma_f^2*exp(-0.5/(l^2)*(x1'*ones(size(x2))...
    - ones(size(x1))'*x2).^2);
%-----\COVARIANCE FUNCTION ----- %
    
K = K_maker(X, X, sigma_f);                 % Training points covariances
K_star = K_maker(x_star, x_star, sigma_f);  % Test point covariances
k_star = K_maker(X, x_star, sigma_f);       % Related covariances

K_full = K + diag(sigma_n^2*ones(1,n_obs));
K_inv = inv(K_full);                        % Inversion

f_star = k_star'*K_inv*y;                   % Outputs

V = K_star - k_star'*K_inv*k_star;          % Full covariance

% log_marginal = -0.5*(y'*K_inv*y - log(det(K_full)) - numel(y)*log(2*pi));

dVdz = diff(f_star)/(x_star(2)-x_star(1));

for i = 1:n_test-1
	if dVdz(i) > beta_required
		min_alt = x_star(i);
		
		for j = i+1:n_test-1
			if dVdz(j) < beta_required
				max_alt = x_star(j);
				break;
			elseif j == n_test-1;
				max_alt = x_star(j);
			end
		end
		if max_alt - min_alt > 3
			break;
		end
	elseif i == n_test-1
		min_alt = 0;
		max_alt = 0;
	end
end


%% Plotting
yellow = [0.6, .8, 1];
red = [1, 0.6, 0.6];

% if any(~(get(0, 'Children')-99))
% 	set(0,'CurrentFigure',99);
% else
	figure(99);
% end
clf; hold on;
sigma_y = sqrt(diag(V));
for i = 1:length(sigma_y)
    if imag(sigma_y(i)) > 1e-6
        fprintf(1, ['\nImaginary component of covariance in element %i',...
            ' with magnitude %5f\n'], i, imag(sigma_y(i)));
    elseif imag(sigma_y(i)) > 0
        fprintf(1, ['\nImaginary component of covariance in element %i',...
            ' with magnitude %1.4e has been ignored\n'], ...
            i, imag(sigma_y(i)));
        sigma_y(i) = real(sigma_y(i));
    else
        
    end
end

fill([f_star'+2*sigma_y', reverse(f_star'-2*sigma_y')], ...
    [x_star, flipdim(x_star, 2)], yellow, 'edgecolor', yellow.^2);
% fill([f_star'+1*sigma_y', reverse(f_star'-1*sigma_y')], ...
%     [x_star, flipdim(x_star, 2)], red, 'edgecolor', red.^2);

plot(f_star', x_star, '-b', 'linewidth', 1.5);
plot(y, X, 'k.')

real_data = wind_field([zeros(2, n_test); -x_star]);
real_data = sqrt(sum(real_data.^2, 1));
plot(real_data, x_star, 'r--', 'linewidth', 1.0);

if min_alt == 0 && max_alt == 0
	plot([y(1), y(end)], [obs_loc(1), obs_loc(1)], 'g--');
else
	plot([y(1), y(end)], [min_alt, min_alt], '--', 'color', [0 .8 0]);
	plot([y(1), y(end)], [max_alt, max_alt], '--', 'color', [0 .8 0]);
end

title('Gaussian Process Fit of Wind Profile');
xlabel('Horizontal Wind Speed, \it{W_x} \rm(m/s)'); ylabel('Altitude (m)');
legend({'2\sigma bounds', 'Fitted curve', 'Data points'...
    , 'True profile', 'Altitude limits'}, 'Location', 'NorthWest');
grid on;


	
	
	