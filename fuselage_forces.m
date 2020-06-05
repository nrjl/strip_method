function [FForces, FMoments] = fuselage_forces(X)
global PLANE_PARAM

Ceb = calc_Ceb(X(7:9));

uvw_wind = X(1:3) - Ceb'*wind_field(X(10:12));
airspeed = sqrt(sum(uvw_wind.^2));

d = PLANE_PARAM.fuse_d;
l = PLANE_PARAM.fuse_l;

[T, P, rho, a, mu, nu] = atmos(-X(12), 0);		% Atmospheric properties
Re = airspeed*l/nu;

Cf = 0.427/((log(Re) - 0.407)^2.64);
Cd = Cf*(1+1.5*(d/l)^3/2 + 7*(d/l)^3) + PLANE_PARAM.Cd0;	% Total drag coefficient

D = Cd*0.5*rho*airspeed^2*(pi*d^2/4);

u_hat = uvw_wind/airspeed;

FForces = -u_hat*D;

r = Ceb*[-l*(1-PLANE_PARAM.fuse_cg)/2; 0; 0];
% FMoments = [0;0;0];
FMoments = cross(r, FForces);