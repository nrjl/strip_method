function [X_new] = euler_stupid(X_old, U, dt, ddt)

X_new = zeros(size(X_old, 1), dt/ddt+1);
X_new(:,1) = X_old;

for i = 1:dt/ddt
	X_new(:,i+1) = X_new(:,i) + ddt*state_rates(X_new(:,i), U);
end

X_new = X_new(:,2:end);