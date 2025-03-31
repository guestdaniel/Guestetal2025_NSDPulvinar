% transform_postthalamus_0pt5mm_to_MNI_and_func1mm.m
% 
% Uses nsd_mapdata to transform Mike Arcaro's posterior thalamus labels to MNI space
%
% Note: all NSD data is LPI arrangement; however, when we transform to MNI the output is RPI arragement.
% Therefore, before saving the outputs to disk with nsd_savenifti we flip them to LPI coordinates.
%
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
combined_map = zeros(182, 218, 182, 8);
for subj=1:8
	% Do MNI first
	combined_map(:, :, :, subj) = nsd_mapdata(...
		subj, 'anat0pt5', 'MNI', ...
		['/home/stone/generic/Dropbox/nsdthalamus/mike/postthalamus/subj0' num2str(subj) '/thalamus.nii.gz'], ...
		'wta', 0, ...
		['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/mni/postthalamus.nii.gz']);
	% Do func1mm next
	nsd_mapdata(...
		subj, 'anat0pt5', 'func1pt0', ...
		['/home/stone/generic/Dropbox/nsdthalamus/mike/postthalamus/subj0' num2str(subj) '/thalamus.nii.gz'], ...
		'wta', 0, ...
		['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/func1mm/postthalamus.nii.gz']);
end
% <<< NOTE >>> nsd_mapdata will return volumes in RPI arrangement because we requested MNI (per Kendrick 04/22/2020)
% Find average ROI across subjects by taking mode and save to disk
combined_map = mode(combined_map, 4);
% Flip combined map to transform it to LPI
combined_map = flip(combined_map, 1);
% Save map to disk
nsd_savenifti(combined_map, [1 1 1], '/home/surly-raid3/dguest-data/Guestetal2021_data/group/mni/postthalamus.nii.gz', 1, [92, 127, 73]);
% <<< NOTE >>> nsd_savenifti assumes that volumes are LPI, but since we're saving in MNI space we should align to MNI template origin of [92, 127, 73]

