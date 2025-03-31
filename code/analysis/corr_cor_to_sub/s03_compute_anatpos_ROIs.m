% analyses/corr_cor_to_sub/04_compute_anatpos_ROIs.m
%
% This script computes the average anatomical position (i.e., A-P, L-M, S-I) for each ROI used in the
% corr_cor_to_sub analysis. These anatomical position data are then saved to disk.

config_Guestetal2025_NSDPulvinar;
hemis = {'lh', 'rh'};

% Load ROI
load([data_dir 'group/fsaverage/compiled_cortical_rois.mat']);

% Load anatpos maps
anatpos = zeros(163842, 3, 2, 8);
names = {'anatpos_anteriorposterior.mgz', 'anatpos_lateralmedial.mgz', 'anatpos_superiorinferior.mgz'};
for subj=1:8
	for hemi=1:2
		for type=1:3
			ap(:, type, hemi, subj) = cvnloadmgz(fullfile(data_dir, ['subj0' num2str(subj)], 'fsaverage', [hemis{hemi} '.' names{type}]));
		end
	end
end

% Compute averages within each ROI
results = zeros(3, 2, 8, 14);
for subj=1:8
	for hemi=1:2
		for type=1:3
			for roi_label=1:14
				results(type, hemi, subj, roi_label) = ...
					mean(ap(compiled_roi(:, hemi, subj, roi_label) == 1, type, hemi, subj));
			end
		end
	end
end

% Verify that results are reasonable visually
% Extract one example map
%ex = zeros(163842, 1);
%for roi_label=1:14
%	ex(compiled_roi(:, 1, 1, roi_label) == 1) = results(1, 1, 1, roi_label);
%end
%cvnlookup('fsaverage', 13, ex);

% Save to disk
save(fullfile(data_dir, 'group', 'cortical_roi_average_positions.mat'), 'results');
