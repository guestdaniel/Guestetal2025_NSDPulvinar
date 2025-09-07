% analyses/prf/cluster_data.m
%
% This script takes the subcortical NSD beta weights and attempts to cluster them 
% using a combination of PCA (for dimensionality reduction) followed by k-means 
% clustering. This is an alternate approach to characterizing the functional responses 
% in the pulvinar suggested by a reviewer.
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];
hemis = {'lh', 'rh'};

% Loop through subjects
for subj=1:8
    for idx_hemi=[1]
        % Print progress
        fprintf('Subj = %d Hemi = %s\n', subj, hemis{idx_hemi})

        % Load subcortical data
        load(['/home/surly-raid3/dguest-data/subcortical/functional/betas/' 'subj' sprintf('%02d', subj) '_subcortical_betas_MNI_fast.mat']);

        % Define LH/RH indices for MNI space (full volume / full brain)
        LH = 1:floor(182/2);
        RH = (max(LH)+1):182;

        % Define LH/RH indices for MNI subspace containing pulvinar data spanning from coords_start(1):coords_end(1), etc...
        LH_sub = LH(LH >= coords_start(1)) - coords_start(1) + 1;
        RH_sub = RH(RH <= coords_end(1)) - coords_start(1) + 1;

        % Reshape data
        x = storage_ind; 
        x = permute(x, [2, 3, 4, 1, 5]);  % convert to (n_x=56, n_y=22, n_z=27, n_sess=n_sess_per_subj(subj), n_trial=750)
        if idx_hemi == 1
            x = x(LH_sub, :, :, :, :);
        else
            x = x(RH_sub, :, :, :, :);
        end
        x = squish(x, 3); % convert to (n_voxel, n_sess, 750)
        x = zscore(x, 0, 3);  % z-score each session
        x = permute(x, [1, 3, 2]); % convert to (n_voxel, 750, n_sess)
        x = reshape(x, [size(x, 1), size(x, 2)*size(x, 3)]); % convert to (n_voxel, n_trial) roughly (~15k, 30k)

        % Use PCA to reduce X from full size (~30k, 30k) to reduced size (~30k, 200)
        % Note that 200 is an arbitrary choice; full decision should involve a scree plot
        [c, s, l] = pca(x, 'NumComponents', 200);   % returns (coefs, score, latent)
        % Note matrix s is size (n_voxel, n_component)

        % Apply k-means clustering to the PCA loadings/scores
        for k = [2, 5, 10, 20]
            cluster_ids = kmeans(s, k, 'MaxIter', 1000, 'Replicates', 10);  % kmeans takes (n_voxel, n_component) -> (n_voxel, cluster_id)

            % Embed data brick in MNI space, branching based on whether we're filling in the
            % LH or RH results based on the current value of `hemi`
            vol = nan(182, 218, 182);	
            if idx_hemi == 1
                subvol = zeros(length(LH_sub), 22, 27);
            else
                subvol = zeros(length(RH_sub), 22, 27);
            end
            subvol(:) = cluster_ids(:);
            if idx_hemi == 1
                vol(LH_sub, ...
                    coords_start(2):coords_end(2), ...
                    coords_start(3):coords_end(3)) = subvol;
            else
                vol(RH_sub, ...
                    coords_start(2):coords_end(2), ...
                    coords_start(3):coords_end(3)) = subvol;
            end

            % Save to disk as nifti 
            fn = fullfile(data_dir, ['subj0' num2str(subj)], 'mni', ['beta_clusters_kmeans_hemi=' hemis{idx_hemi} '_k=' num2str(k) '.nii.gz']);
            nsd_savenifti(vol, [1, 1, 1], fn, 1, [92, 127, 73]);
        end
    end
end