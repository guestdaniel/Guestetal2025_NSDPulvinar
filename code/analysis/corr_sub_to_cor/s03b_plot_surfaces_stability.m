% analyses/corr_sub_to_cor/s02b_plot_surfaces_stability.m
%
% Plots the results of `s02b_compute_correlation_stability.m`
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Set up a few things
seeds = {'contrast_Nmax_lh', 'contrast_Nmax_rh', 'body_Nmax_lh', 'body_Nmax_rh', 'contrast_Nmin_lh', 'contrast_Nmin_rh', 'body_Nmin_lh', 'body_Nmin_rh'};
cutoffs = {0.05*j, 0.025*j, 0.02*j, 0.01*j, []};
N = 20;  % how many top N voxels we analyze

% Load in average visual cortex rois
lh_roi = cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'lh.prf-visualrois_probavg.mgz'));
rh_roi = cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'rh.prf-visualrois_probavg.mgz'));
roi_visual_cortex = [lh_roi; rh_roi];	

% Load in floc rois
roi_bodies = [cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'lh.floc-bodies_probavg.mgz'));
	      cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'rh.floc-bodies_probavg.mgz'))];
roi_faces = [cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'lh.floc-faces_probavg.mgz'));
	      cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'rh.floc-faces_probavg.mgz'))];
roi_places = [cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'lh.floc-places_probavg.mgz'));
	      cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'rh.floc-places_probavg.mgz'))];
roi_words = [cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'lh.floc-words_probavg.mgz'));
	      cvnloadmgz(fullfile(data_dir, 'group', 'fsaverage', 'rh.floc-words_probavg.mgz'))];

% % Step 1: Handle the standard maps (single correlation using all trials per
% % vertex). We plot a few different map flavors (e.g., mean, t-value
% % thresholded, etc.)
% % Loop through methods
for method=1:2
	% Loop through seeds
	for idx_seed=1:8
        % Loop through subjects and load all correlation maps
        for n=1:N
            % Determine paths to data
            data = {};
            % Loop through subjects and load results
            for subj=1:8
                fn = fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', ['corr_sub_to_cor_' seeds{idx_seed} '_voxel' num2str(n) '_method_' num2str(method) '.mat']);
                load(fn);
                data{subj} = corr_matrix;
            end

            % Average across data
            data_mean = nanmean(cell2mat(data), 2);

            % Plot maps
            % Here, to avoid excessive indenting, we put several for loops at the same indentation depth, each explained below
            for cutoff=1:length(cutoffs)   % Loop through thresholding variants (either at 0.01 or no thresholding)
            for label=1:3 % Loop through whether to plot visual cortex roi labels
                % Handle label
                if label == 1
                    options = {'hemibordercolor', 'w', 'rgbnan', 1, ...
                        'roimask', {roi_visual_cortex == 1,
                                roi_visual_cortex == 2,
                                roi_visual_cortex == 3,
                                roi_visual_cortex == 4,
                                roi_visual_cortex == 5,
                                roi_visual_cortex == 6,
                                roi_visual_cortex == 7}, ...
                        'roicolor', [255, 255, 255]/255, ...
                        'roiwidth', 1};
                elseif label == 2
                    options = {'hemibordercolor', 'w', 'rgbnan', 1, ...
                        'roimask', {roi_bodies,
                                roi_faces,
                                roi_words,
                                roi_places}, ...
                        'roicolor', {[0, 0, 255]/255, 
                                [0, 255, 0]/255, 
                                [0, 255, 255]/255,
                                [255, 255, 255]/255}, ...
                        'roiwidth', 1};
                else
                    options = {'hemibordercolor', NaN, 'rgbnan', 1};
                end

                % Set clims
                clims = [-0.15, 0.15];

                % Plot surface plot of means
                [rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_mean, clims, cmapsign4(256), cutoffs{cutoff}, [], 0, options);
                if cutoff == length(cutoffs)
                    th = 'none';
                else
                    th = num2str(100*imag(cutoffs{cutoff}));
                end
                fn = fullfile(plot_dir_bulk, ['sub_to_cor_map_supplemental_seed=' seeds{idx_seed} '_method=' num2str(method) '_threshold=' th '_label=' num2str(label) '_voxel=' num2str(n) '.png']);
                imwrite(rgbimg, fn);

                % Display
                disp('Running!')
            end
            end
        end
	end
end