% analyses/corr_sub_to_cor/01_compute_correlations.m
%
% Computes correlations between responses in the "contrast peak" and "person peak" voxels in the subcortex and responses over the whole cortical surface 
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Parameters
n_sess_per_subj = [40, 40, 32, 30, 40, 32, 40, 30];
skip_mode = 1; % if 1 and file exists, don't repeat analysis for that subj/cond

% Identify target voxels for each subject (body peak and contrast peak)
[cont_lh, cont_rh, body_lh, body_rh] = f_identify_candidate_voxels_subject_native();

% Labels for the four possible seed ,ocations
labels = {'contrast_lh', 'contrast_rh', 'body_lh', 'body_rh'};

% Loop through subjects and compute correlations
for subj=1:8
	% Load data
	[data_cortical, data_subcortical] = f_load_data_subject_native(subj, n_sess_per_subj(subj), 1); 

	% Load this subject's mat file indicating where data in brick belongs in subject-native volume
	load(fullfile(xxx_kkdata_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));
	coords_start = [d1(1), d2(1), d3(1)];
	coords_end = [d1(end), d2(end), d3(end)];
	%save([output_dir '/subj0' num2str(subj) '_coords_start.mat'], 'coords_start');
	%save([output_dir '/subj0' num2str(subj) '_coords_end.mat'], 'coords_end');
	data_subcortical = reshape(data_subcortical, length(d1), length(d2), length(d3), size(data_subcortical, 2));

	% Create array of coordinates and array of labels
	coords = {cont_lh(subj, :), cont_rh(subj, :), body_lh(subj, :), body_rh(subj, :)};
	labels = {'contrast_lh', 'contrast_rh', 'body_lh', 'body_rh'};

	% Loop over coords and labels and compute each correlation a single
	% time (for baseline analysis... bootstrapped variant is provided
	% below!)
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
		corr_matrix = nanmean(corrs, 2);
		% Save result to disk
		save(fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_2.mat']), 'corr_matrix');
    end
end

% Before we run analysis below, we need to fetch 
d = {};
for subj = 1:8
	% Load this subject's mat file indicating where data in brick belongs in subject-native volume
	load(fullfile(xxx_kkdata_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));
	%coords_start = [d1(1), d2(1), d3(1)];
	%coords_end = [d1(end), d2(end), d3(end)];
    d{subj} = {d1, d2, d3};
end

out_method1 = cell(8, 1);
out_method2 = cell(8, 1);
 
parpool(4);
for subj = 1:8  
	% Loop over coords and labels and compute each correlation with
	% boostrapping. To do so, we more or less follow the logic as above
    % but immediately before each calculation we randomly sample with
    % replacement new row/obs indices and use those instead of true row/obs
    % indices. This process is repeated many times, and then we save a
    % matrix of each correlation value based on each repeat of the
    % algorithm. Thus, the output correlation matrix sizes are (n_vertex,
    % n_rep_boot) instead of (n_vertex, 1).
    
    % Fetch task so we can print task ID
    % t = getCurrentTask();
    
    % Set params
    n_rep_boot = 1000;  % how many times to repeat bootstrap
    
    % Load data
	[data_cortical, data_subcortical] = f_load_data_subject_native(subj, n_sess_per_subj(subj), 1);
    d1 = d{subj}{1};
    d2 = d{subj}{2};
    d3 = d{subj}{3};
	data_subcortical = reshape(data_subcortical, length(d1), length(d2), length(d3), size(data_subcortical, 2));

	% Create array of coordinates and array of labels
	coords = {cont_lh(subj, :), cont_rh(subj, :), body_lh(subj, :), body_rh(subj, :)};
    
    % Create empty cell arrays to store results
    out_method1{subj} = cell(2, 1);
    out_method2{subj} = cell(2, 1);
    
    % Loop over seeds
	for ii=[2 4]
		% Extract coords we're using
		c_this = coords{ii};

		% ===== METHOD #1 =====
        % Determine file name
        fn = fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_1_bootstrap.mat']);
        
        % Skip this analysis if file exists
        if (isfile(fn) && (skip_mode == 1))
            % Print status
            fprintf("Skipping subj %d, seed %d, method 1 since it already exists!\n", subj, ii);
        else
            % Print status
            fprintf("Running subj %d, seed %d, method 1 since it does not exist!\n", subj, ii);
            
            % Pre-allocate storage for results
            n_vertex = size(data_cortical, 1);
            n_trial = size(data_cortical, 2);
            corr_matrix = zeros(n_vertex, n_rep_boot);

            % Subset data to only what we strictly need in subcortex
            temp_sub = squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), :))';

            % Loop over reps and calculate correlation, sampling rows with
            % replacement
            tic;
            for rr=1:n_rep_boot
                % Display progress
                if mod(rr, 5) == 0
                    fprintf("Bootstrapping Subj %d, Seed %d, Method 1... Step %d\n", subj, ii, rr);
                end
                trialidxs = randi(n_trial, n_trial, 1);
                corr_matrix(:, rr) = f_calc_corr(temp_sub(trialidxs), data_cortical(:, trialidxs))';        
            end

            out_method1{subj}{ii/2} = corr_matrix;
        end
        
		% ===== METHOD #2 =====
        % Determine filename
        fn = fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_2_bootstrap.mat']);
        
        % Skip this analysis if file exists
        if (isfile(fn) && (skip_mode == 1))
            % Print status
            fprintf("Skipping subj %d, seed %d, method 2 since it already exists!\n", subj, ii);
        else
            % Print status
            fprintf("Running subj %d, seed %d, method 2 since it does not exist!\n", subj, ii);
            
            % Identify possible pairings of subcortical and cortical indices
            [trial_idxs_subcortical, trial_idxs_cortical] = f_select_trials_method2(subj, n_sess_per_subj(subj));

            % Create empty matrix for correlations 
            corrs = zeros(size(corr_matrix, 1), n_rep_boot, 6);

            % Subset subcortical down to only what we need
            temp_sub = squeeze(data_subcortical(c_this(1), c_this(2), c_this(3), :))';

            % Loop over reps and calculate correlation, resampling rows/obs
            % with replacement
            for rr=1:n_rep_boot
                % Display progress
                if mod(rr, 5) == 0
                    fprintf("Bootstrapping Subj %d, Seed %d, Method 2... Step %d\n", subj, ii, rr);
                end

                % Draw row/obs indices and use to resample trial_idxs_... from
                % above for this step of the bootstrap
                trialidxs = randi(size(trial_idxs_subcortical, 1), size(trial_idxs_subcortical, 1), 1);
                tx_sub_resamp = trial_idxs_subcortical(trialidxs, :);
                tx_cor_resamp = trial_idxs_cortical(trialidxs, :);

                % Loop over sub/cor trial configurations
                for repeat=1:6
                    idxs_sub = tx_sub_resamp(:, repeat);
                    idxs_cor = tx_cor_resamp(:, repeat);
                    idxs_sub = idxs_sub(~isnan(idxs_sub));
                    idxs_cor = idxs_cor(~isnan(idxs_cor));
                    corrs(:, rr, repeat) = f_calc_corr(temp_sub(idxs_sub), data_cortical(:, idxs_cor))';
                end
            end
            % Average across repeats
            corr_matrix = nanmean(corrs, 3);
            
            % Save result to disk
            out_method2{subj}{ii/2} = corr_matrix;
        end
	end
end

% Loop to save
% Loop over subjs
for subj = 1:8
    % Only subloop if not empty
    if ~isempty(out_method1{subj})
        % Loop over seed locations
        for ii = [2 4]
            % Method 1
            fn = fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_1_bootstrap.mat']);
            if ~isempty(out_method1{subj}{ii/2})
                corr_matrix = out_method1{subj}{ii/2};
                save(fn, 'corr_matrix');
            end
        end
    end
    
    % Only subloop for method2 if data structure is not empty
    if ~isempty(out_method2{subj})
        % Loop over seed locations
        for ii = [2 4]
            % Method 2
            fn = fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' labels{ii} '_method_2_bootstrap.mat']);
            if ~isempty(out_method2{subj}{ii/2})
                corr_matrix = out_method2{subj}{ii/2};
                save(fn, 'corr_matrix');
            end
        end
    end
end