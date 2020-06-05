aa = dir('DS_runs')

for ii = 1:numel(aa)
	bb = dir(['DS_runs\', aa(ii).name, '\*.mat']);
	for jj = 1:numel(bb)
		if ~isempty(strfind(bb(jj).name, 'dynamic'))
			load(['DS_runs\', aa(ii).name, '\', bb(jj).name]);
			if exist('modetimes', 'var')
				fprintf(1, '\n %s     ', ['DS_runs\', aa(ii).name, '\', bb(jj).name]);
% 				fprintf(1, '%0.3f ', (modetimes(1,:)));
				modetimes
			end
		end
	end
end