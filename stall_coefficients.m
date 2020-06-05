
% phi_le = 0;
% r_le = 0;
% phi_te = 0;
alfa = alfa_full*pi/180;

Cd_90 = 1.7 + (0.3 - phi_le*(0.2+0.08*phi_le))*(1-1.8*sqrt(r_le))...
	- phi_te*(0.2+0.08*phi_te);

Cn = Cd_90*sin(alfa)./(0.56 + 0.44*sin(alfa));

gamma = 0.28*sqrt(r_le);

Ct = 0.5*0.0075*cos(alfa) + Cn.*sin(alfa);

ClStC = 1.5*(Cn.*cos(alfa) - Ct.*sin(alfa));
CdStC = Cn.*sin(alfa) + Ct.*cos(alfa);