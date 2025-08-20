% analyses/prf/cluster_data.m
%
% This script takes the subcortical NSD beta weights and attempts to cluster them 
% using a combination of PCA (for dimensionality reduction) followed by k-means 
% clustering.
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];

% Loop through subjects
for subj=1:8
	% Load subcortical data
	load(['/home/surly-raid3/dguest-data/subcortical/functional/betas/' 'subj' sprintf('%02d', subj) '_subcortical_betas_MNI_fast.mat']);
	x = storage_ind;
	x = permute(x, [2, 3, 4, 1, 5]);  % convert to (n_x, n_y, n_z, n_sess, 750)
	x = squish(x, 3); % convert to (n_voxel, n_sess, 750)
	x = zscore(x, 0, 3);  % z-score each session
	x = permute(x, [1, 3, 2]); % convert to (n_voxel, 750, n_sess)
	x = reshape(x, [size(x, 1), size(x, 2)*size(x, 3)]); % convert to (n_voxel, n_trial) roughly (~30k, 30k)

    % Use PCA to reduce X from full size (~30k, 30k) to reduced size (~30k, 200)
    % Note that 200 is an arbitrary choice; full decision should involve a scree plot
    [c, s, l] = pca(x, 'NumComponents', 200);   % returns (coefs, score, latent)
    % Note matrix s is size (n_voxel, n_component)

    % Apply k-means clustering to the PCA loadings/scores
    for k = [2, 5, 10, 20, 50]
        cluster_ids = kmeans(s, k);  % kmeans takes (n_voxel, n_component) -> (n_voxel, cluster_id)

        % Embed data brick in MNI space
        vol = zeros(182, 218, 182);	
		subvol = zeros(56, 22, 27);
        subvol(:) = cluster_ids(:);
        vol(coords_start(1):coords_end(1), ...
            coords_start(2):coords_end(2), ...
            coords_start(3):coords_end(3)) = subvol;

        % Save to disk as nifti 
        fn = fullfile(data_dir, ['subj0' num2str(subj)], 'mni', ['beta_clusters_kmeans_k=' num2str(k) '.nii.gz']);
        nsd_savenifti(vol, [1, 1, 1], fn, 1, [92, 127, 73]);
    end
end