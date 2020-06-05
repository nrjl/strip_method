load s2048_control.mat
s2048_control = {polars, Re, def};

alfa = -90:1:90;
Re = 150000;
def = -5:1:5;

Cl = zeros(length(def), length(alfa));
Cd = Cl;
Cm = Cl;

legen = cell(1, length(def));

for i = 1:length(def)	
	for j = 1:length(alfa)
		
		[Cl(i,j), Cd(i,j), Cm(i,j)] = ...
			wing_lookup(alfa(j), Re, def(i), s2048_control);
		
	end
	legen{i} = sprintf('def = %0.3g', def(i));
end

%%
figure(9); clf;
plot(alfa, Cl);
xlabel('Angle of attack (deg)'); ylabel('Cl');
legend(legen);

figure(10); clf;
plot(alfa, Cd);
xlabel('Angle of attack (deg)'); ylabel('Cd');
legend(legen);

figure(11); clf;
plot(alfa, Cm);
xlabel('Angle of attack (deg)'); ylabel('Cm');
legend(legen);
