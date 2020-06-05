function K = squared_exp(x1,x2, sigma_f, l)
%----- COVARIANCE FUNCTION ----- %
K = sigma_f^2*exp(-0.5/(l^2)*(x1'*ones(size(x2)) - ones(size(x1))'*x2).^2);
%-----\COVARIANCE FUNCTION ----- %