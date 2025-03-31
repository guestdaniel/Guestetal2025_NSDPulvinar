% transform_thalamus_anat0pt8_to_MNI.m
% 
% Uses nsd_mapdata to transform hand-drawn thalamus ROI labels 0.8mm isotropic anatomical space to MNI space.
%
% Note: all NSD data is LPI arrangement; however, when we transform to MNI the output is RPI arragement.
% Therefore, before saving the outputs to disk with nsd_savenifti we flip them to LPI coordinates.
%
% Use nsd_mapdata to map each subject's subcortical ROIs to new ROIs in MNI space
combined_map = zeros(182, 218, 182, 8);
for subj=1:8
	combined_map(:, :, :, subj) = nsd_mapdata(...
		subj, 'anat0pt8', 'MNI', ...
		['/home/surly-raid4/kendrick-data/nsd/nsddata/ppdata/subj0' num2str(subj) '/anat/roi/thalamus.nii.gz'], ...
		'wta', 0, ...
		['/home/surly-raid3/dguest-data/Guestetal2021_data/subj0' num2str(subj) '/mni/thalamus.nii.gz']);
end
% <<< NOTE >>> nsd_mapdata will return volumes in RPI arrangement because we requested MNI (per Kendrick 04/22/2020)
% Find average ROI across subjects by taking mode and save to disk
combined_map = mode(combined_map, 4);
% Flip combined map to transform it to LPI
combined_map = flip(combined_map, 1);
% Save map to disk
nsd_savenifti(combined_map, [1 1 1], '/home/surly-raid3/dguest-data/Guestetal2021_data/group/mni/thalamus.nii.gz', 1, [92, 127, 73]);
% <<< NOTE >>> nsd_savenifti assumes that volumes are LPI, but since we're saving in MNI space we should align to MNI template origin of [92, 127, 73]
