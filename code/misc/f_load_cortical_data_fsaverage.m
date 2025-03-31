function data_cortical = f_load_cortical_data_fsaverage(subj, n_sess, hemisphere)
% Loads preprocessed subject BOLD data, whole cortex fsaverage 
% 
% Args:
% 	subj (int): which subject to load
%	n_sess (int): how many sessions to load
%	hemisphere (int): which hemisphere to load, 1 = LH, 2 = RH
%
% Returns:
% 	data (array): array of z-scored betas of shape (n_vertex, n_sample)

% Run config and set up directories
config_Guestetal2025_NSDPulvinar
data_directory = [nsdbeta_dir '/ppdata/subj0' num2str(subj) '/fsaverage/betas_fithrf_GLMdenoise_RR/'];

% Load cortical data
data_cortical = zeros(163842, n_sess*750);
if hemisphere == 1
	hemi = 'lh';
else
	hemi = 'rh';
end
for session=1:n_sess
	disp(['Loading cortical data for subj ' num2str(subj) ', session ' num2str(session)]);
	try 
		data = cvnloadmgz([data_directory '/' hemi '.betas_session' sprintf('%02d', session) '.mgh']);
		data = squeeze(data);
		data = zscore(data, 0, 2);  % z-score each session
		data_cortical(:, 750*(session-1)+(1:750)) = data;
	catch ME
		disp(['No session #' num2str(session) ' for subject #' num2str(subj)]);
		continue;
	end
end
