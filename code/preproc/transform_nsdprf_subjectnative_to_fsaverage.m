% transform_prf_func1mm_R2_to_fsvarage.m
% 
% Uses nsd_mapdata to transform pRF R2 output from subject-native functional space to fsaverage space

% Set source/target directory
store_dir = '/home/surly-raid3/dguest-data/Guestetal2021_data/';
nsdfaceprf_dir = '/home/stone-ext4/generic/DropboxOBSOLETE/nsdfaceprf/freesurfer/';

% Select models to transform
models = {'backgroundauto', 'bodyauto', 'contrastNEW', 'faceauto', 'foregroundauto', 'salience', 'wordauto'};

% Use nsd_mapdata to map each subject's pRF results to MNI space
for model=1:length(models)
	for subj=1:8
        % First, load surface in subj-native space
	    temp = cvnloadmgz(fullfile(nsdfaceprf_dir, ['subj0' num2str(subj)], 'label', ['lh.' models{model} '_R2.mgz']));

        % Next, transform to fsaverage space
		out = nsd_mapdata(subj, 'lh.white', 'fsaverage', temp);

        % Finally, save to disk
		nsd_savemgz(out, [store_dir 'subj0' num2str(subj) '/fsaverage/' 'lh.' models{model} '_R2.mgz'], ['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/fsaverage']);

        % Print success
        fprintf('Transformed %s for subject %d to fsaverage space.\n', models{model}, subj);
	end
end

