function [contrast_maximum_lh, contrast_maximum_rh, body_maximum_lh, body_maximum_rh] = f_identify_candidate_voxels_subject_native()
% Locates peaks in the subcortical variance explained for contrast model and body model and returns 
% their spatial locations
%
% Returns:
%	contrast_maximum_lh (array): spatial location of contrast peak in left hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, 3)
%	contrast_maximum_rh (array): spatial location of contrast peak in right hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, 3)
%	body_maximum_lh(array): spatial location of bodyauto peak in left hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, 3)
%	body_maximum_rh(array): spatial location of bodyauto peak in right hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, 3)
config_Guestetal2025_NSDPulvinar;

% Create storage for final outpus
contrast_maximum_lh = zeros(8, 3);
contrast_maximum_rh = zeros(8, 3);
body_maximum_lh = zeros(8, 3);
body_maximum_rh = zeros(8, 3);

% Loop through subjects and load each map
for subj=1:8 
	% Load in posterior thalamus ROI
	roi = load_untouch_nii(fullfile(data_dir, ['subj0' num2str(subj)], 'func1mm', 'postthalamus.nii.gz'));
	roi = roi.img;

	% Load in this subject's .mat file (which contains incides d1, d2, d3
	% indicating where data in brick belongs in the subject-native volume)
	load(fullfile(addtl_dir, ['datab3nativesurface_subj0' num2str(subj) '.mat']));

	% Load in Contrast R2
	temp = load_untouch_nii(fullfile(data_dir, ['subj0' num2str(subj)], 'func1mm', 'contrastNEW_R2.nii.gz'));
	temp = temp.img;
	temp = squeeze(temp(:, :, :, 1));
	% Make empty volume the same size as roi
	vol = zeros(size(roi));
	% Insert data into volume
	vol(d1, d2, d3) = temp;
	% Zero out voxels not in ROI
	vol(roi ~= 1) = 0;
	R2_contrast = vol;

	% Load in body R2
	temp = load_untouch_nii(fullfile(data_dir, ['subj0' num2str(subj)], 'func1mm', 'bodyauto_R2.nii.gz'));
	temp = temp.img;
	temp = squeeze(temp(:, :, :, 1));
	% Make empty volume the same size as roi
	vol = zeros(size(roi));
	% Insert data into volume
	vol(d1, d2, d3) = temp;
	% Zero out voxels not in ROI
	vol(roi ~= 1) = 0;
	R2_bodyauto = vol;

	% Identify LH and RH indices
	LH = 1:floor(size(R2_contrast, 1)/2);
	RH = (max(LH)+1):size(R2_contrast);
	xx_LH = length(LH);
	xx_RH = length(RH);
	yy = size(R2_contrast, 2);
	zz = size(R2_contrast, 3);

	% Locate contrast LH maximum
	temp = R2_contrast;
	temp(RH, :, :) = 0;
	[~, argmax] = max(temp(:));
	[x, y, z] = ind2sub(size(temp), argmax);
	contrast_maximum_lh(subj, :) = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];

	% Locate contrast RH maximum
	temp = R2_contrast;
	temp(LH, :, :) = 0;
	[~, argmax] = max(temp(:));
	[x, y, z] = ind2sub(size(temp), argmax);
	contrast_maximum_rh(subj, :) = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];

	% Locate body LH maximum
	temp = R2_bodyauto;
	temp(RH, :, :) = 0;
	[~, argmax] = max(temp(:));
	[x, y, z] = ind2sub(size(temp), argmax);
	body_maximum_lh(subj, :) = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];

	% Locate body RH maximum
	temp = R2_bodyauto;
	temp(LH, :, :) = 0;
	[~, argmax] = max(temp(:));
	[x, y, z] = ind2sub(size(temp), argmax);
	body_maximum_rh(subj, :) = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];
end

end
