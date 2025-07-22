function [contrast_maximum_lh, contrast_maximum_rh, body_maximum_lh, body_maximum_rh, contrast_minimum_lh, contrast_minimum_rh, body_minimum_lh, body_minimum_rh] = f_identify_candidate_voxels_subject_native_all(N)
% Locates peaks in the subcortical variance explained for contrast model and body model and returns 
% their spatial locations; in contrast to the simpler `f_identify_candidate_voxels_subject_native`,
% which returns only a peak voxel location, this function returns the locations of the top N and
% bottom N voxels.
%
% Returns:
%	contrast_maximum_lh (array): spatial location of contrast top N voxels in left hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, N, 3)
%	contrast_maximum_rh (array): spatial location of contrast top N voxels in right hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, N, 3)
%	body_maximum_lh(array): spatial location of bodyauto top N voxels in left hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, N, 3)
%	body_maximum_rh(array): spatial location of bodyauto top N voxels in right hemisphere, in units 
%		of voxels in subcortical brick containing subcortical ROIs (1mm isotropic resolution), 
%		shape of (n_subj, N, 3)
%   .... same for contrast_minimum_* and body_minimum_* except for bottom N voxels
config_Guestetal2025_NSDPulvinar;

% Create storage for final outputs
contrast_maximum_lh = zeros(8, N, 3);
contrast_maximum_rh = zeros(8, N, 3);
body_maximum_lh = zeros(8, N, 3);
body_maximum_rh = zeros(8, N, 3);
contrast_minimum_lh = zeros(8, N, 3);
contrast_minimum_rh = zeros(8, N, 3);
body_minimum_lh = zeros(8, N, 3);
body_minimum_rh = zeros(8, N, 3);

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

	% Locate best N voxels
    bestN = cell(2, 2);

    % Loop over data sources [contrast, bodyauto]
    for idx_data = 1:2
        % Loop over hemispheres [left, right]
        for idx_hemi = 1:2
            % Based on idx_data, decide if we are processing contrast or body data
            if idx_data == 1
                temp = R2_contrast;
            else
                temp = R2_bodyauto;
            end
            % Based on idx_hemi, decide if we are zeroing out LH or RH
            if idx_hemi == 1  % we're looking at LH, ignoring RH
                temp(RH, :, :) = 0;
            else              % we're looking at RH, ignoring LH
                temp(LH, :, :) = 0;
            end

            % Identify top N voxels and their locations
            [~, idx_sort] = sort(temp(:));
            idx_sort = flipud(idx_sort);  % voxels now sorted with best at idx_sort(1) and worst at idx_sort(end)
            [x, y, z] = ind2sub(size(temp), idx_sort(1:N));
            bestN{idx_data, idx_hemi} = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];
        end
    end

	% Locate worst N voxels
    worstN = cell(2, 2);

    % For worst N, we put every value outside the ROI to a large number 
	R2_contrast(roi ~= 1) = 1000;
	R2_bodyauto(roi ~= 1) = 1000;

    % Loop over data sources [contrast, bodyauto]
    for idx_data = 1:2
        % Loop over hemispheres [left, right]
        for idx_hemi = 1:2
            % Based on idx_data, decide if we are processing contrast or body data
            if idx_data == 1
                temp = R2_contrast;
            else
                temp = R2_bodyauto;
            end
            % Based on idx_hemi, decide if we are filling with big values LH or RH
            if idx_hemi == 1  % we're looking at LH, ignoring RH
                temp(RH, :, :) = 1000;
            else              % we're looking at RH, ignoring LH
                temp(LH, :, :) = 1000;
            end

            % Identify top N voxels and their locations
            [~, idx_sort] = sort(temp(:));
            [x, y, z] = ind2sub(size(temp), idx_sort(1:N));
            worstN{idx_data, idx_hemi} = [x-d1(1)+1, y-d2(1)+1, z-d3(1)+1];
        end
    end

    % Store results for this subject
    contrast_maximum_lh(subj, :, :) = bestN{1, 1};
    contrast_maximum_rh(subj, :, :) = bestN{1, 2};
    body_maximum_lh(subj, :, :) =     bestN{2, 1};
    body_maximum_rh(subj, :, :) =     bestN{2, 2};
    contrast_minimum_lh(subj, :, :) = worstN{1, 1};
    contrast_minimum_rh(subj, :, :) = worstN{1, 2};
    body_minimum_lh(subj, :, :) =     worstN{2, 1};
    body_minimum_rh(subj, :, :) =     worstN{2, 2};
end

end
