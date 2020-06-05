function [AForces, AMoments, varargout] = strip_method(X, U, varargin)
%--------------------------------------------------------------------------
%
% FUNCTION:		strip_method
%
% PURPOSE:		calculate aircraft body forces and moments using a simple
%				strip method
%               
% SYNTAX:		[AForces, AMoments] = strip_method(X, U)
%				[AForces, AMoments, a_handle] = strip_method(X, U, 1)
%				The second syntax will plot the force vector arrows and
%				return a handle to them in the third output argument
%
% INPUTS:		X		- state vector
%				U		- control vector
%
% OUTPUTS:		AForces	- Aerodynamic forces
%				AMoments- Aerodynamic moments
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		October 2007
%
% MODIFIED:     October 2007
%
% See also:		strip_forces, state_rates
%--------------------------------------------------------------------------
global PLANE_AERO

if nargin == 3
	plot_flag = varargin{1};
else
	plot_flag = 0;
end

% stime = cputime;

n_surf = length(PLANE_AERO);
l_conv = eye(4);
l_conv(2,2) = -1;

n_poa = 0;
for i = 1:n_surf
	n_poa = n_poa + ...
		(size(PLANE_AERO(i).qc_surf, 2)-1)*(1 + PLANE_AERO(i).mirror);
end

poa = zeros(3, n_poa);
F = poa;
M = poa;
stalled = zeros(1, n_poa);
cn = 1;


for i = 1:n_surf
	qc_surf = PLANE_AERO(i).qc_surf;	% Quarter chord and chord lengths
	base_alfa = PLANE_AERO(i).base_alfa(:);	% Base angle of incidence
	mirror = PLANE_AERO(i).mirror;		% Mirrored surface about yz plane
	profile = i;						% Surface number coefficients
	control = PLANE_AERO(i).control(:);	% Control number
	reverse = PLANE_AERO(i).reverse(:);	% Controls reversed
	c_lim = PLANE_AERO(i).c_lim; % * [1, 0; 0, -1];
	
	U2 = [0; U];					% Augmented control vector so that index 1 becomes uncontrolled
	c_lim2 = [[0;0], c_lim];		% Augmented control limits (same as above)
	
	% Deflections
	reverse_vector = -2*(reverse<0)+1;	% Index -1 for control reversal, 1 for no reversal
	control_vector = U2(control.*(control>0)+1).*(reverse_vector);
	def = (control_vector > 0).*control_vector.*c_lim2(1,abs(control)+1)' ...
		- (control_vector < 0).*control_vector.*c_lim2(2,abs(control)+1)';
	
	% Angle of incidence variations
	control_vector = U2(-control.*(control<0)+1).*(reverse_vector);
	ai = base_alfa + (control_vector > 0).*control_vector.*c_lim2(1,abs(control)+1)' ...
		- (control_vector < 0).*control_vector.*c_lim2(2,abs(control)+1)';
	
	fn = cn + size(PLANE_AERO(i).qc_surf, 2)-2;
	
	% Calculate strip forces
	[poa(:, cn:fn), F(:, cn:fn), M(:, cn:fn), stalled(cn:fn)] = ...
		strip_forces(qc_surf, X, ai, @wing_lookup, profile, def, mirror);
	
	cn = fn+1;
			
	if mirror
		qc_surf = l_conv*qc_surf;
		reverse_vector = ones(size(reverse_vector)) - 2*(reverse > 0);

		control_vector = U2(control.*(control>0)+1).*(reverse_vector);
		def = (control_vector > 0).*control_vector.*c_lim2(1,abs(control)+1)' ...
			- (control_vector < 0).*control_vector.*c_lim2(2,abs(control)+1)';

		% Angle of incidence variations
		control_vector = U2(-control.*(control<0)+1).*(reverse_vector);
		ai = base_alfa + (control_vector > 0).*control_vector.*c_lim2(1,abs(control)+1)' ...
			- (control_vector < 0).*control_vector.*c_lim2(2,abs(control)+1)';

		fn = cn + size(PLANE_AERO(i).qc_surf, 2)-2;

		% Calculate strip forces
		[poa(:, cn:fn), F(:, cn:fn), M(:, cn:fn), stalled(cn:fn)] = ...
			strip_forces(qc_surf, X, ai, @wing_lookup, profile, def, 1);
		
		cn = fn+1;	
	end
	
end

if plot_flag
	F_scale = 0.1;
	Ceb = calc_Ceb(X(7:9));

	F_all = Ceb*F*F_scale;
	F_end_nstall = poa + F_all.*([1;1;1]*~stalled);
	F_end_stall = poa + F_all.*([1;1;1]*stalled);

	c_ax = axis;
	c_dis = sqrt(sum(c_ax([2,4,6]).^2 - c_ax([1,3,5]).^2));
	
	ah = sqrt(sum(F_all.^2, 1))*0.15*72/c_dis/3;	% Arrow height (25%)
	if any(~stalled)
		a_handle{1} = arrow3(poa', (F_end_nstall)', 'b-', ah/3, ah, max(ah)/8);
	else
		a_handle{1} = 0;
	end
	
	if any(stalled)
		a_handle{2} = arrow3(poa', (F_end_stall)', 'r-', ah/3, ah, max(ah)/8);
	else
		a_handle{2} = 0;
	end
	varargout(1) = {a_handle};
end

AForces = sum(F, 2);
AForces = AForces.*~(abs(AForces) < 1e-10);

AMoments = sum(M, 2);
AMoments = AMoments.*~(abs(AMoments) < 1e-9);