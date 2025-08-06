% analyses/corr_cor_to_sub/s03b_compute_correlations_partial.m
%
% Computes correlations between responses in a block containing the visual thalamus and surrounding tissue and responses in average cortical seeds.
% This version of the script is based on the first round of reviews and attempts to use partial correlation to isolate unique contributions of cortical seed
% areas.
%
% Depends
%	group/fsaverage/compiled_cortical_betas.mat
%
% Outputs:
%	subj*/mni/corr_cor_to_sub_*.nii.gz

% Load cortical data and set number of sessions per subject
config_Guestetal2025_NSDPulvinar;
load([data_dir 'group/fsaverage/compiled_cortical_betas.mat']);  % size (30000, 2, 8, 14) for (n_trial, n_hemi, n_subj, n_seed)
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];

% Loop through subjects, hemispheres, and cortical ROI labels
for subj=1:8
	% Load subcortical data and reshape to desired format
	load(['/home/surly-raid3/dguest-data/subcortical/functional/betas/' 'subj' sprintf('%02d', subj) '_subcortical_betas_MNI_fast.mat']);
	x = storage_ind;
	x = permute(x, [2, 3, 4, 1, 5]);  % convert to (n_x, n_y, n_z, n_sess, 750)
	x = squish(x, 3); % convert to (n_voxel, n_sess, 750)
	x = zscore(x, 0, 3);  % z-score each session
	x = permute(x, [1, 3, 2]); % convert to (n_voxel, 750, n_sess)
	x = reshape(x, [size(x, 1), size(x, 2)*size(x, 3)]); % conver to (n_voxel, n_trial)

    % Loop over hemi
	for hemi=1:2
        % First, subset the cortical data to only be this subject and hemisphere, giving (30000, 14) matrix
        data_cortical = squeeze(avg_data(:, hemi, subj, :));

        % Next, loop over label
		for label=1:14
            % Display current state
			disp(['Subj: ' num2str(subj) ', Hemi: ' num2str(hemi) ', Label: ' num2str(label)]);

            % Break cortical data into seed of interest (indexed by label) and other seeds, which will be treated as regressors in partial correlation
            data_cortical_target = squeeze(data_cortical(:, label));
            data_cortical_regressors = data_cortical(:, 1:14 ~= label);

			% ===== METHOD #1 =====
			% Calculate partial correlations coefficients
                % subcortical input is (n_trial, n_voxel)
                % cortical target input is (n_trial, 1)
                % cortical regressors input is (n_trial, 13)
			corr = f_calc_partial_corr(x(:, 1:(750*n_sess_per_subj(subj)))', data_cortical_target(1:(750*n_sess_per_subj(subj))), data_cortical_regressors(1:(750*n_sess_per_subj(subj)), :));  

            % Create empty volume and subvolume to store results
			subvol = zeros(56, 22, 27);
			vol = zeros(182, 218, 182);	

			% Embed data brick in MNI space
			subvol(:) = corr(:);
			vol(coords_start(1):coords_end(1), ...
			    coords_start(2):coords_end(2), ...
			    coords_start(3):coords_end(3)) = subvol;

			% Save to disk as nifti 
			nsd_savenifti(vol, [1, 1, 1], [data_dir 'subj0' num2str(subj) '/mni/corr_cor_to_sub_hemi_' num2str(hemi) '_label_' num2str(label) '_method_1_partial.nii.gz'], 1, [92, 127, 73]);

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
				corrs(:, repeat) = ...
                    f_calc_partial_corr(x(:, idxs_sub)', data_cortical_target(idxs_cor), data_cortical_regressors(idxs_cor, :));  
			end
			corr = nanmean(corrs, 2);

			% Create empty volume and subvolume to store results
			subvol = zeros(56, 22, 27);
			vol = zeros(182, 218, 182);	

			% Embed data brick in MNI space
			subvol(:) = corr(:);
			vol(coords_start(1):coords_end(1), ...
			    coords_start(2):coords_end(2), ...
			    coords_start(3):coords_end(3)) = subvol;

			% Save to disk as nifti 
			nsd_savenifti(vol, [1, 1, 1], [data_dir 'subj0' num2str(subj) '/mni/corr_cor_to_sub_hemi_' num2str(hemi) '_label_' num2str(label) '_method_2_partial.nii.gz'], 1, [92, 127, 73]);
		end
	end
end
