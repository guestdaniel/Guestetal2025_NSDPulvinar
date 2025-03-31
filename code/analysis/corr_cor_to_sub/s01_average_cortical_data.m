% Config
config_Guestetal2025_NSDPulvinar;

% Define number of sessions for each subject
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];

% Load compiled rois
load([data_dir 'group/fsaverage/compiled_cortical_rois.mat']); % Compiled ROI is size (n_vertex, n_hemi, n_subj, n_label)

% Storage variables
avg_data = zeros(30000, 2, 8, 14);  % (30000, n_hemi, n_subj, n_label)

% Loop through subjects and hemispheres
for subj=1:8
	for hemi=1:2
		% Load data from disk (already z-score
		data = load_cortical_data_fsaverage(subj, n_sess_per_subj(subj), hemi);
		% Average data within ROIs
		for label=1:14
			avg_data(1:(750*n_sess_per_subj(subj)), hemi, subj, label) = nanmean(data(compiled_roi(:, hemi, subj, label) == 1, :), 1);
		end
	end
end

% Save to disk
save([data_dir 'group/fsaverage/compiled_cortical_betas.mat'], 'avg_data');
