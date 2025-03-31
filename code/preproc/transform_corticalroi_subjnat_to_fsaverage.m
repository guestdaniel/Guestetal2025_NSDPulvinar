% transform_cortical_roi_subjnat_to_fsaverage.m
% 
% Uses nsd_mapdata to transform visual cortex ROIs and high-level visual ROIs to fsaverage from 
% subject-native surface space
%
% ===============================================
% = Map visual cortex ROIs to fsaverage         =
% ===============================================
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
for subj=1:8
	hemis = {'lh', 'rh'};
	for hemi=1:2
		output = nsd_mapdata(...  
			subj, [hemis{hemi} '.white'], 'fsaverage', ...
			['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj0' num2str(subj) '/label/' hemis{hemi} '.prf-visualrois.mgz'], ...
			'wta', 0);
		nsd_savemgz(output, ['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/fsaverage/' hemis{hemi} '._prf-visualrois.mgz'], ...
		['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/fsaverage']);
	end
end

% ===============================================
% = Map high-level cortex ROIs to fsaverage     =
% ===============================================
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
rois = {'flocbodiestval', 'flocfacestval', 'flocplacestval', 'flocwordtval', 'floc-bodies', 'floc-faces', 'floc-places', 'floc-words'};
for subj=1:8
	hemis = {'lh', 'rh'};
	for hemi=1:2
		for roi=1:8
			output = nsd_mapdata(...
				subj, [hemis{hemi} '.white'], 'fsaverage', ...
				['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj0' ...
				 num2str(subj) '/label/' hemis{hemi} '.' rois{roi} '.mgz'], ...
				'wta', 0);
			nsd_savemgz(output, ['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.' rois{roi} '.mgz'], ...
				    ['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/fsaverage']);
         	    end
	end
end

% ===============================================
% = Map Kastner/Wang atlas ROIs to fsaverage    =
% ===============================================
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
for subj=1:8
	hemis = {'lh', 'rh'};
	for hemi=1:2
		output = nsd_mapdata(...
			subj, [hemis{hemi} '.white'], 'fsaverage', ...
			['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/subj0' ...
			 num2str(subj) '/label/' hemis{hemi} '.' 'Kastner2015' '.mgz'], ...
			'wta', 0);
		nsd_savemgz(output, ['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.' 'Kastner2015' '.mgz'], ...
			    ['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/fsaverage']);
	end
end
