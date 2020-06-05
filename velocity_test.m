%% Velocity test
R = 2;
ri = 2;
vc = 2.5;

r = 1;

t = linspace(-pi, 0, 50);

vz = -vc * (1 - (R + r*cos(t))/(R + ri)).*cos(t);

figure(1); clf; hold on;
plot(R + r*cos(t), vz);
plot(R + r*cos(t(25)), vz(25), 'rx')


%% Linear piecewise
x = linspace(R - 3*ri, R + 3*ri, 1001);
dex = find((x > R-ri) .* (x < R+ri));

vz2 = (x <= R).* (vc - x*vc/R) + ...
      (x >  R).* (-vc*(x - R).*(2*R - x))./x/R;

%% Cubic
vz3 = vc/(2*R^2)*(x.^2 - 3*R*x + 2*R^2);

%% Sinusoid
vz4 = R*vc/pi*sin(pi/R*x);

%% Quintic
% k = 2.8;
% c1 = (vc - 2*k)/(2*R^4);
% c2 = (4*k - vc)/(2*R^4);
% c3 = -k;

call = rref([R^5,    R^3,   R, 0;...
             5*R^4,  3*R^2, 1, vc;...
             10*R^2, 3,     0, 0]);
c1 = call(1,4); c2 = call(2,4); c3 = call(3,4);

vz5 = c1*(x-R).^5 + c2*(x-R).^3 + c3*(x-R);

%% Seventh order
call = rref([R^6,      R^4,     R^2,    1, 0 ;...
             7*R^6,    5*R^4,   3*R^2,  1, vc;...
             21*R^4,   10*R^3,  3,      0, 0 ;...
             1216*R^6, 208*R^4, 28*R^2, 1, 0 ]);
c1 = call(1,5); c2 = call(2,5); c3 = call(3,5); c4 = call(4,5); 

vz6 = c1*(x-R).^7 + c2*(x-R).^5 + c3*(x-R).^3 + c4*(x-R);

%% Plotting
figure(2); clf;

subplot(2,2,1); hold on;
plot([0, 1, 2], [0, 0, 0], 'r.', x/R, vz2.*x, 'b-')    
plot(x/R, vz3.*x, 'k--')
plot(x/R, vz4, 'g-.')
plot(x/R, vz5, 'r:')
axis([min(x/R), max(x/R), -15, 15]);
% plot(x, vz5, 'm-')

legend({'required points', 'linear piecewise', 'quadratic', ...
    'sinusoid', 'quartic'})
title('\it{h(x-R)}')

subplot(2,2,2); hold on;
plot([0, 1, 2], [vc, 0, 0], 'r.', x/R, vz2, 'b-')    
plot(x/R, vz3, 'k--')
plot(x/R, vz4./x, 'g-.')
plot(x/R, vz5./x, 'r:')
% plot(x, vz6./x, 'm-')
axis([min(x/R), max(x/R), -5, 5]);
title('\it{f(x)}')
xlabel('\it{}x/R');
ylabel('\it{}V_z/V_c');

subplot(2,2,3); hold on;
plot([0, 1, 2], [0, 0, 0], 'r.', x(dex)/R, vz2(dex).*x(dex), 'b-')
plot(x(dex)/R, vz3(dex).*x(dex), 'k--')
plot(x(dex)/R, vz4(dex), 'g-.')
plot(x(dex)/R, vz5(dex), 'r:')
% plot(x(dex), vz6(dex), 'm-')
% axis([R - ri, R + ri, -vc/2, vc/2])
title('\it{h(x-R)} \rm area of interest')

subplot(2,2,4); hold on;
plot([0, 1, 2], [vc, 0, 0], 'r.', x(dex)/R, vz2(dex), 'b-')    
plot(x(dex)/R, vz3(dex), 'k--')
plot(x(dex)/R, vz4(dex)./x(dex), 'g-.')
plot(x(dex)/R, vz5(dex)./x(dex), 'r:')
% plot(x(dex), vz6(dex)./x(dex), 'm-')
% axis([R - ri, R + ri, -vc/2, vc])
title('\it{f(x)} \rm area of interest')
xlabel('\it{}x/R');
ylabel('\it{}V_z/V_c');