% transform_anatpos_subjnat_to_fsaverage.m
% 
% Uses nsd_mapdata to transform anatomical position maps from subject-native surface space to fsaverage
%
% ===============================================
% = Map visual cortex ROIs to fsaverage         =
% ===============================================
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
names = {'anatpos_anteriorposterior.mgz', 'anatpos_lateralmedial.mgz', 'anatpos_superiorinferior.mgz'};
for subj=1:8
	hemis = {'lh', 'rh'};
	for hemi=1:2
		for type=1:3
			% First, load the surface
			x = cvnloadmgz(['/home/stone/generic/Dropbox/nsdfaceprf/freesurfer/subj0' num2str(subj) '/label/' hemis{hemi} '.APRLSI.mgz']);
			% Next, map the surface to fsaverage
			output = nsd_mapdata(subj, [hemis{hemi} '.white'], 'fsaverage', x(:, 1, 1, type), 'nearest', 0);
			% Finally, save output to disk!
			nsd_savemgz(output, ['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.' names{type}], ['/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/fsaverage']);
		end
	end
end

