% analyses/corr_cor_to_sub/03_compute_correlations.m
%
% Computes correlations between responses in a block containing the visual thalamus and surrounding tissue and responses in average cortical seeds
%
% Depends
%	group/fsaverage/compiled_cortical_betas.mat
%
% Outputs:
%	subj*/mni/corr_cor_to_sub_*.nii.gz

% Load cortical data and set number of sessions per subject
config_Guestetal2025_NSDPulvinar;
load([data_dir 'group/fsaverage/compiled_cortical_betas.mat']);
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];

% Loop through subjects, hemispheres, and cortical ROI labels
for subj=1:8
	% Load subcortical data
	load(['/home/surly-raid3/dguest-data/subcortical/functional/betas/' 'subj' sprintf('%02d', subj) '_subcortical_betas_MNI_fast.mat']);
	x = storage_ind;
	x = permute(x, [2, 3, 4, 1, 5]);  % convert to (n_x, n_y, n_z, n_sess, 750)
	x = squish(x, 3); % convert to (n_voxel, n_sess, 750)
	x = zscore(x, 0, 3);  % z-score each session
	x = permute(x, [1, 3, 2]); % convert to (n_voxel, 750, n_sess)
	x = reshape(x, [size(x, 1), size(x, 2)*size(x, 3)]); % conver to (n_voxel, n_trial)
	for hemi=1:2
		for label=1:14
			disp(['Subj: ' num2str(subj) ', Hemi: ' num2str(hemi) ', Label: ' num2str(label)]);
			% ===== METHOD #1 =====
			% Calculate correlations
			corr = f_calc_corr(x(:, 1:(750*n_sess_per_subj(subj))), ...
					 squeeze(avg_data(1:(750*n_sess_per_subj(subj)), hemi, subj, label))');
			% create empty MNI volume
			subvol = zeros(56, 22, 27);
			vol = zeros(182, 218, 182);	

			% Embed data brick in MNI space
			subvol(:) = corr(:);
			vol(coords_start(1):coords_end(1), ...
			    coords_start(2):coords_end(2), ...
			    coords_start(3):coords_end(3)) = subvol;

			% Save to disk as nifti 
			nsd_savenifti(vol, [1, 1, 1], [data_dir 'subj0' num2str(subj) '/mni/corr_cor_to_sub_hemi_' num2str(hemi) '_label_' num2str(label) '_method_1.nii.gz'], 1, [92, 127, 73]);

			% ===== METHOD #2 =====
			[trial_idxs_subcortical, trial_idxs_cortical] = f_select_trials_method2(subj, n_sess_per_subj(subj));

			% Create empty matrix for correlations
			corrs = zeros(size(x, 1), 6); 
			% Calculate correlations
			for repeat=1:6
				idxs_sub = trial_idxs_subcortical(:, repeat);
				idxs_cor = trial_idxs_cortical(:, repeat);
				idxs_sub = idxs_sub(~isnan(idxs_sub));
				idxs_cor = idxs_cor(~isnan(idxs_cor));
				corrs(:, repeat) = f_calc_corr(x(:, idxs_sub), ...
						 squeeze(avg_data(idxs_cor, hemi, subj, label))');
			end
			corr = nanmean(corrs, 2);

			% create empty MNI volume
			subvol = zeros(56, 22, 27);
			vol = zeros(182, 218, 182);	

			% Embed data brick in MNI space
			subvol(:) = corr(:);
			vol(coords_start(1):coords_end(1), ...
			    coords_start(2):coords_end(2), ...
			    coords_start(3):coords_end(3)) = subvol;

			% Save to disk as nifti 
			nsd_savenifti(vol, [1, 1, 1], [data_dir 'subj0' num2str(subj) '/mni/corr_cor_to_sub_hemi_' num2str(hemi) '_label_' num2str(label) '_method_2.nii.gz'], 1, [92, 127, 73]);
		end
	end
end
