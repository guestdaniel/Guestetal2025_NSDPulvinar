% analyses/corr_sub_to_cor/02_plot_surfaces.m
%
% Plots surface data derived from 01_compute_correlations.m
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;

% Set up a few things
labels = {'contrast_lh', 'contrast_rh', 'body_lh', 'body_rh'};
variants = {0.01*j, []};

% Load in average visual cortex rois
lh_roi = cvnloadmgz(['~/subcortical/rois/lh.prf-visualrois_probavg.mgz']);
rh_roi = cvnloadmgz(['~/subcortical/rois/rh.prf-visualrois_probavg.mgz']);
roi_visual_cortex = [lh_roi; rh_roi];	
% Load in floc rois
roi_bodies = {};
roi_faces = {};
roi_objects = {};
roi_words = {};
roi_bodies = [cvnloadmgz(['~/subcortical/rois/lh.floc-bodies_probavg.mgz']);
	      cvnloadmgz(['~/subcortical/rois/rh.floc-bodies_probavg.mgz'])];
roi_faces = [cvnloadmgz(['~/subcortical/rois/lh.floc-faces_probavg.mgz']);
	      cvnloadmgz(['~/subcortical/rois/rh.floc-faces_probavg.mgz'])];
roi_places = [cvnloadmgz(['~/subcortical/rois/lh.floc-places_probavg.mgz']);
	      cvnloadmgz(['~/subcortical/rois/rh.floc-places_probavg.mgz'])];
roi_words = [cvnloadmgz(['~/subcortical/rois/lh.floc-words_probavg.mgz']);
	      cvnloadmgz(['~/subcortical/rois/rh.floc-words_probavg.mgz'])];

% Loop through methods and labels
for method=1:2
	for label=1:4
		data = {};
		% Loop through subjects and load results
		for subj=1:8
			load([data_dir 'subj0' num2str(subj) '/fsaverage/' 'corr_sub_to_cor_' labels{label} '_method_' num2str(method) '.mat']);
			data{subj} = corr_matrix;
		end
		% Average across data
		data = nanmean(cell2mat(data), 2);
		for var=1:1   % Loop through thresholding variants (either at 0.01 or no thresholding)
		for label=1:3 % Loop through whether to plot visual cortex roi labels
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
				options = {'hemibordercolor', 'w', 'rgbnan', 1};
			end
			% Plot surface plot
			[rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data, ...
				[-0.15 0.15], cmapsign4(256), variants{var}, [], 0, options);
			imshow(rgbimg);
			saveas(gcf, fullfile(plot_dir_bulk, [labels{label} '_' num2str(method) '_' num2str(var) '_' num2str(label) '.png']));
			close;
		end
		end
		end

	end
end

