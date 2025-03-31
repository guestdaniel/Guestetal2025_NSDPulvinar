% analyses/corr_sub_to_cor/01_compute_correlations.m
%
% Computes correlations between responses in the "contrast peak" and "person peak" voxels in the subcortex and responses over the whole cortical surface 
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Parameters
n_voxel_along_axis = 10;
%n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];
n_sess_per_subj = [3, 3, 3, 3, 3, 3, 3, 3];

% Identify target voxels for each subject (body peak and contrast peak)
[cont_lh, cont_rh, body_lh, body_rh] = f_identify_candidate_voxels_subject_native();

% Loop through subjects and compute correlations
for subj=1:8
	% Load data
	[data_cortical, data_subcortical] = f_load_data_subject_native(subj, n_sess_per_subj(subj), 1); 

	% Load this subject's mat file indicating where data in brick belongs in subject-native volume
	load(fullfile(addtl_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));
	coords_start = [d1(1), d2(1), d3(1)];
	coords_end = [d1(end), d2(end), d3(end)];
	%save([output_dir '/subj0' num2str(subj) '_coords_start.mat'], 'coords_start');
	%save([output_dir '/subj0' num2str(subj) '_coords_end.mat'], 'coords_end');
	data_subcortical = reshape(data_subcortical, length(d1), length(d2), length(d3), size(data_subcortical, 2));

	% Create array of coordinates and array of labels
	coords = {cont_lh(subj, :), cont_rh(subj, :), body_lh(subj, :), body_rh(subj, :)};
	labels = {'contrast_lh', 'contrast_rh', 'body_lh', 'body_rh'};

	% Loop over coords and labels
	for ii=1:4
		% Extract coords we're using
		c_this = coords{ii};

		% ===== METHOD #1 =====
		% Compute correlation (left hemisphere subcortical, both hemispheres cortical)
		corr_matrix = f_calc_corr(squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), :))', ...
			   		  data_cortical)';
		% Save result to disk
		save(fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_1.mat']), 'corr_matrix');
		% ===== METHOD #2 =====
		[trial_idxs_subcortical, trial_idxs_cortical] = f_select_trials_method2(subj, n_sess_per_subj(subj));
		% Create empty matrix for correlations 
		corrs = zeros(size(corr_matrix, 1), 6);
		% Loop through repeats and perform analysis
		for repeat=1:6
			idxs_sub = trial_idxs_subcortical(:, repeat);
			idxs_cor = trial_idxs_cortical(:, repeat);
			idxs_sub = idxs_sub(~isnan(idxs_sub));
			idxs_cor = idxs_cor(~isnan(idxs_cor));
			corrs(:, repeat) = f_calc_corr(squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), idxs_sub))', ...
					 	       data_cortical(:, idxs_cor))';

		end
		% Average across repeats
		corr = nanmean(corrs, 2);
		% Save result to disk
		save([data_dir 'subj0' num2str(subj) '/fsaverage/' 'corr_sub_to_cor_' labels{ii} '_method_2.mat'], 'corr_matrix');

	end
end
