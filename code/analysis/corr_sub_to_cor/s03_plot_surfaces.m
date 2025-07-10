% analyses/corr_sub_to_cor/02_plot_surfaces.m
%
% Plots surface data derived from 01_compute_correlations.m
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Set up a few things
seeds = {'contrast_lh', 'contrast_rh', 'body_lh', 'body_rh'};
variants = {0.05*j, 0.025*j, 0.02*j, 0.01*j, []};
t_criterion = 3.499;  % for df=7, this is the two-sided t-val for p=0.01

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
% for method=1:2
% 	% Loop through seeds
% 	for seed=[2 4]
% 		data = {};
% 		% Loop through subjects and load results
% 		for subj=1:8
% 			load([data_dir 'subj0' num2str(subj) '/fsaverage/' 'corr_sub_to_cor_' seeds{seed} '_method_' num2str(method) '.mat']);
% 			data{subj} = corr_matrix;
% 		end
% 		% Average across data
% 		data_mean = nanmean(cell2mat(data), 2);
% 		data_std =  std(cell2mat(data), 0, 2);
% 		data_t = data_mean ./ (data_std/sqrt(8));
% 
% 		% Visualize
% 		for var=1:length(variants)   % Loop through thresholding variants (either at 0.01 or no thresholding)
% 		for label=1:3 % Loop through whether to plot visual cortex roi labels
% 		for crange=1:2 % Loop through whether we pick a limited color range or one normalized to 99th percentile
% 			% Handle label
% 			if label == 1
% 				options = {'hemibordercolor', 'w', 'rgbnan', 1, ...
% 					   'roimask', {roi_visual_cortex == 1,
% 						       roi_visual_cortex == 2,
% 						       roi_visual_cortex == 3,
% 						       roi_visual_cortex == 4,
% 						       roi_visual_cortex == 5,
% 						       roi_visual_cortex == 6,
% 						       roi_visual_cortex == 7}, ...
% 					   'roicolor', [255, 255, 255]/255, ...
% 					   'roiwidth', 1};
% 			elseif label == 2
% 				options = {'hemibordercolor', 'w', 'rgbnan', 1, ...
% 					   'roimask', {roi_bodies,
% 						       roi_faces,
% 						       roi_words,
% 						       roi_places}, ...
% 					   'roicolor', {[0, 0, 255]/255, 
% 							[0, 255, 0]/255, 
% 							[0, 255, 255]/255,
% 							[255, 255, 255]/255}, ...
% 					   'roiwidth', 1};
% 			else
% 				options = {'hemibordercolor', NaN, 'rgbnan', 1};
% 			end
% 
% 			% Handle crange
% 			if crange == 1
% 				clims = [-0.15, 0.15];
% 			else
% 				percentile = quantile(data_mean, 0.99);
% 				clims = [-percentile, percentile];
% 			end
% 
% 			% Plot surface plot of means
% 			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_mean, ...
% 				clims, cmapsign4(256), variants{var}, [], 0, options);
% 			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' num2str(var) '_' num2str(label) '_' num2str(crange) '.png']));
% 
% 			% Plot surface plot of means (shuffled in space)
% 			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_mean(randperm(length(data_mean))), ...
% 				clims, cmapsign4(256), variants{var}, [], 0, options);
% 			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' num2str(var) '_' num2str(label) '_' num2str(crange) '_shuffled.png']));
% 
% 			% Plot surface plot of t
% 			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_t, ...
% 				[-10 10], cmapsign4(256), t_criterion*j, [], 0, options);
% 			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' num2str(var) '_' num2str(label) '_' num2str(crange) '_tval.png']));
% 
% 			% Plot surface plot of means, thresholded by t-val instead of arbitrary cutoff
% 			temp = data_mean;
% 			temp(abs(data_t) < t_criterion) = 0;
% 			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, temp, ...
% 				clims, cmapsign4(256), 0.0*j, [], 0, options);
% 			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' num2str(var) '_' num2str(label) '_' num2str(crange) '_threshold_by_tval.png']));
% 
% 			% Display
% 			disp('Running!')
% 		end
% 		end
%         end
% 	end
% end

% Step 2: Handle the boostrapped maps!
% For the boostrapped maps, we don't try various thresholds for
% visualization; we either plot it unthresholded, or we plot it thresholded
% based on a bootstrapped significance test looking to see that at least
% 95% of the replicates have an r value that exceeds 0.

% Loop through methods
for method=1:2
	% Loop through seeds
	for seed=[2 4]
		data = [];
		% Loop through subjects and load results
		for subj=1:8
            if method == 1
                load([data_dir 'subj0' num2str(subj) '/fsaverage/' 'corr_sub_to_cor_' seeds{seed} '_method_' num2str(method) '_bootstrap.mat']);
            else
                load([data_dir 'subj0' num2str(subj) '/fsaverage/' 'corr_sub_to_cor_' seeds{seed} '_method_' num2str(method) '_bootstrap.mat']);
            end
            if subj == 1
                data = corr_matrix;
            else
                data = cat(3, data, corr_matrix);
            end
		end
		% Average across data
        n_rep_boot = size(data, 2);
		data_mean = nanmean(nanmean(data, 3), 2);
        sig_by_subj = squeeze((sum(data > 0.0, 2)/n_rep_boot) > 0.95);
        sig_net = sum(sig_by_subj, 2) >= 5;

		% Visualize
		for label=1:3 % Loop through whether to plot visual cortex roi labels
		for crange=1:2 % Loop through whether we pick a limited color range or one normalized to 99th percentile
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

			% Handle crange
			if crange == 1
				clims = [-0.15, 0.15];
			else
				percentile = quantile(data_mean, 0.99);
				clims = [-percentile, percentile];
            end
           
			% Plot surface plot of means
			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_mean, ...
				clims, cmapsign4(256), [], [], 0, options);
			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' 'na' '_' num2str(label) '_' num2str(crange) '_bootstrapped_no_thr.png']));

            % Plot surface plot of means
            data_filt = data_mean;
            data_filt(sig_net == 0) = 0.0;
			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data_filt, ...
				clims, cmapsign4(256), [0.001], [], 0, options);
			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' 'na' '_' num2str(label) '_' num2str(crange) '_bootstrapped_thr_95_7.png']));
            
            % Plot surface plot of means
			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, sig_net, ...
				clims, cmapsign4(256), [], [], 0, options);
			imwrite(rgbimg, fullfile(plot_dir_bulk, [seeds{seed} '_' num2str(method) '_' 'na' '_' num2str(label) '_' num2str(crange) '_incl_after_bootstrap.png']));

            
			% Display
			disp('Running!')

		end
        end
	end
end
