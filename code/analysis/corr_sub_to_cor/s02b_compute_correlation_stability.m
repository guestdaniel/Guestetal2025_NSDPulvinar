% s02b_compute_correlation_stability.m
%
% Computes average correlation between particular subcortical voxels' BOLD and
% each vertex in the surface. In contrast to the simpler `s02_compute_correlations.m`,
% this function uses `f_identify_candidate_voxels_subject_native_all.m` to identify
% the best and worst 20 voxels, instead of only the single best voxel, and correlate
% all of them with the cortical surface.
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Parameters
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];
N = 20;  % how many voxels to use in analysis

% Identify target voxels for each subject (body peak and contrast peak)
[cont_lh, cont_rh, body_lh, body_rh, cont_lh_m, cont_rh_m, body_lh_m, body_rh_m] = f_identify_candidate_voxels_subject_native_all(N);

% Labels for the four possible seed locations
labels = {'contrast_Nmax_lh', 'contrast_Nmax_rh', 'body_Nmax_lh', 'body_Nmax_rh', 'contrast_Nmin_lh', 'contrast_Nmin_rh', 'body_Nmin_lh', 'body_Nmin_rh'};

% Loop through subjects and compute correlations
for subj=1:8
	% Load data
	[data_cortical, data_subcortical] = f_load_data_subject_native(subj, n_sess_per_subj(subj), 1); 

	% Load this subject's mat file indicating where data in brick belongs in subject-native volume
	load(fullfile(xxx_kkdata_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));
	coords_start = [d1(1), d2(1), d3(1)];
	coords_end = [d1(end), d2(end), d3(end)];
	data_subcortical = reshape(data_subcortical, length(d1), length(d2), length(d3), size(data_subcortical, 2));

	% Create array of coordinates specifically for this subject
	coords = {cont_lh(subj, :, :), cont_rh(subj, :, :), body_lh(subj, :, :), body_rh(subj, :, :), cont_lh_m(subj, :, :), cont_rh_m(subj, :, :), body_lh_m(subj, :, :), body_rh_m(subj, :, :)};

	% Loop over coords and labels and compute each correlation 
	for ii=1:length(labels) 
        % Extract coordinates for this feature/side's N voxels
        Ncoords = squeeze(coords{ii});

        % Loop over voxels
        for jj=1:N
            % Set C_this to be the 3D coordinate for the particular seed voxel we're examining
            c_this = Ncoords(jj, :);

            % ===== METHOD #1 =====
            % Compute correlation (one hemisphere subcortical, both hemispheres cortical)
            corr_matrix = f_calc_corr(squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), :))', data_cortical)';
            % Save result to disk
            save(fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_voxel' num2str(jj) '_method_1.mat']), 'corr_matrix');
            
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
                corrs(:, repeat) = f_calc_corr(squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), idxs_sub))', data_cortical(:, idxs_cor))';
            end
            % Average across repeats
            corr_matrix = nanmean(corrs, 2);
            % Save result to disk
            save(fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_voxel' num2str(jj) '_method_2.mat']), 'corr_matrix');
        end
    end
end