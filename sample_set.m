function [new_z, new_wind] = sample_set(old_z, old_wind, obs_z, obs_wind, n_samples)
%--------------------------------------------------------------------------
%
% FUNCTION:		sample_set
%
% PURPOSE:		determine whether a new sample should be added to the
%				sample set
%               
% SYNTAX:		[new_z, new_wind] = sample_set(old_z, old_wind, obs_z,
%					obs_wind, n_samples)
%
% INPUTS:		old_z		- existing set of observation locations
%				old_wind	- existing set of wind observations
%				obs_z		- new observation location
%				obs_wind	- new wind observation
%				n_samples	- size of sample set
%
% OUTPUTS:		new_z		- new set of observation locations
%				new_wind	- new set of wind observations
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		July 2008
%
% MODIFIED:      2008
%
% See also:		find_ds_limits
%--------------------------------------------------------------------------
nz = size(old_z, 2);
min_z = min(old_z);
max_z = max(old_z);

if nz < n_samples
	if nz == 0
		new_z = obs_z;
		new_wind = obs_wind;
	else
		temp = sortrows([old_z, obs_z; old_wind, obs_wind]');
		new_z = temp(:,1)';
		new_wind = temp(:,2:end)';
	end
	return
	
elseif obs_z < min_z
	dex = find(old_z == min_z);
	old_z(dex) = obs_z;
	old_wind(:,dex) = obs_wind;
	new_z = old_z;
	new_wind = old_wind;
	return
	
elseif obs_z > max_z
	dex = find(old_z == max_z);
	old_z(dex) = obs_z;
	old_wind(:,dex) = obs_wind;
	new_z = old_z;
	new_wind = old_wind;
	return

else
	zspan = max_z-min_z;
	dex = floor((obs_z - min_z)/zspan*nz)+1;
	tvalue = min_z + (dex-1)/(nz-1)*zspan; %min_z + dex/nz*zspan;
	if abs(obs_z - tvalue) < abs(old_z(dex) - tvalue)
		old_z(dex) = obs_z;
		old_wind(:,dex) = obs_wind;
	end
	
	new_z = old_z;
	new_wind = old_wind;
end