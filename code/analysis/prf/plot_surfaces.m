% analyses/prf/plot_surfaces.m
%
% Plots average variance explained data for surface preps of NSD PRF models
%
% Handle paths and settings 
config_Guestetal2025_NSDPulvinar;
store_dir = '/home/surly-raid3/dguest-data/Guestetal2021_data/';  % TODO: replace with stable paths

% Set up a few things
models = {'backgroundauto', 'bodyauto', 'contrastNEW', 'faceauto', 'foregroundauto', 'salience', 'wordauto'};
datatypes = {'R2'};

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

% Loop through labeling schemes
for model=1:length(models)
    % Loop through subjects and load data for each, then average across subjects
    data = {};
    for subj=1:8
        data{subj} = cvnloadmgz(fullfile(store_dir, ['subj0' num2str(subj)], '/fsaverage/', ['lh.' models{model} '_R2.mgz']));
    end

    % Average across subjects
    if size(cell2mat(data), 4) == 1
        data = nanmean(cell2mat(data), 2);
    else
        temp = cell2mat(data);
        data = nanmean(squeeze(temp(:, :, :, 1)), 2);  % WTF is going on here?
    end

    % Loop over labeling schemes
    for label=1:3
        % Set up labeling and other parameters for options
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
        [rawimg, Lookup, rgbimg, himg] = cvnlookup('fsaverage', 10, data, [0.0 40.0], hot(256), [], [], 0, options);

        % Show and save result
        imshow(rgbimg);
        saveas(gcf, fullfile(plot_dir_bulk, ['fsaverage_' models{model} '_R2_label' num2str(label) '.png']));
        close;

        % Print result
        fprintf('Plotted model %s with label scheme %d.\n', models{model}, label);
    end
end

