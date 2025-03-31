function [data_cortical, data_subcortical] = f_load_data_subject_native(subj, n_sess, hemisphere)
% Loads preprocessed subject BOLD data, whole cortex fsaverage + subject native thalamus
% 
% Args:
% 	subj (int): which subject to load
%	n_sess (int): how many sessions to load
%	hemisphere (int): which hemisphere to load, 1 = LH, 2 = RH
%
% Returns:
% 	data_cortical (array): array of z-scored betas of shape (n_vertex, n_sample)
%	data_subcortical (array): array of z-scored betas of shape (n_voxel, n_sample)

% Set up directories
config_Guestetal2025_NSDPulvinar;

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
		data = cvnloadmgz(fullfile(nsdbeta_dir, 'ppdata', ['subj0' num2str(subj)], 'fsaverage', 'betas_fithrf_GLMdenoise_RR', [hemi '.betas_session' sprintf('%02d', session) '.mgh']));
		data = squeeze(data);
		data = zscore(data, 0, 2);  % z-score each session
		data_cortical(:, 750*(session-1)+(1:750)) = data;
	catch ME
		disp(['No session #' num2str(session) ' for subject #' num2str(subj)]);
		continue;
	end
end

% Load subcortical data
% Load in this subject's .mat file (which contains incides d1, d2, d3
% indicating where data in brick belongs in the subject-native volume)
load(fullfile(addtl_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));

% Load in this subject's func1mm T1 for reference (LPI)
T1 = load_untouch_nii(fullfile(nsd_dir, 'ppdata', ['subj0' num2str(subj)], 'func1mm', 'T1_to_func1mm.nii.gz'));
T1 = T1.img;

% Create empty array of size T1
mask = zeros(size(T1));
mask(min(d1):max(d1), min(d2):max(d2), min(d3):max(d3)) = 1;

% Load subject's functional data within thalamus ROIs in subject native space @ 1mm
disp(['Loading subcortical data for subj ' num2str(subj)]);
data = f_load_voxels(subj, n_sess, 3, 2, 1, 1, mask, 0); % shape: n_sess x n_voxel x n_trial 
% Z-score each session of subcortical functional data 
data = zscore(data, 0, 3);	
% Reshape data data to concatenated in time
data = permute(data, [2 3 1]); % size: n_voxel x n_trial x n_sess
n_trial = size(data, 2);
data_subcortical = reshape(data, size(data, 1), size(data, 2)*size(data, 3)); % size: n_voxel x n_trial*n_sess
end


