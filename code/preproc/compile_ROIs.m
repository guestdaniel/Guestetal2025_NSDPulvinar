% analyses/corr_cor_to_sub/01_compile_ROIs.m
% 
% This script compiles the cortical ROIs needed for the cortical -> subcortical correlation analysis
% These ROIs are:
%	[1] V1
%	[2] V2
%	[3] V3
%	[4] hV4
%	[5] OFA
%	[6] FFA
%	[7] aTL-faces
%	[8] EBA
%	[9] FBA
%  	[10] PPA
%	[11] VWFA
%	[12] MT
%	[13] IPS
%	[14] FEF
% 
% Depends:
%	cortical ROIs in fsaverage (visualrois and floc) --- preproc/transform_corticalroi_subjnat_to_fsaverage.m	
% Outputs:
%	compiled_cortical_rois.mat		

% Configure
config_Guestetal2025_NSDPulvinar;
hemis = {'lh', 'rh'};

% Storage variables
compiled_roi = zeros(163842, 2, 8, 14);  % [n_vertex, n_hemi, n_subj, n_label]

% Loop through subjects and hemispheres
for subj=1:8
	for hemi=1:2
		% Visual areas
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '._' 'prf-visualrois' '.mgz']);
		compiled_roi(:, hemi, subj, 1) = (roi == 1 | roi == 2);  % V1
		compiled_roi(:, hemi, subj, 2) = (roi == 3 | roi == 4);  % V2
		compiled_roi(:, hemi, subj, 3) = (roi == 5 | roi == 6);  % V3
		compiled_roi(:, hemi, subj, 4) = (roi == 7);             % hV4
		% Face areas
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.floc-faces' '.mgz']);
		compiled_roi(:, hemi, subj, 5) = (roi == 1);             % OFA
		compiled_roi(:, hemi, subj, 6) = (roi == 2 | roi == 3);  % FFA
		compiled_roi(:, hemi, subj, 7) = (roi == 5);             % aTL-faces
		% Body areas 
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.floc-bodies' '.mgz']);
		compiled_roi(:, hemi, subj, 8) = (roi == 1);             % EBA
		compiled_roi(:, hemi, subj, 9) = (roi == 2 | roi == 3);  % FBA
		% Place areas
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.floc-places' '.mgz']);
		compiled_roi(:, hemi, subj, 10) = (roi == 2);             % PPA
		% Word areas
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.floc-words' '.mgz']);
		compiled_roi(:, hemi, subj, 11) = (roi == 2 | roi == 3);  % VWFA
		% Kastner ROIs 
		roi = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.Kastner2015' '.mgz']);
		compiled_roi(:, hemi, subj, 12) = (roi == 13 | roi == 12);  % TO1/TO2 (area MT)
		compiled_roi(:, hemi, subj, 13) = (roi >= 19 & roi <= 23);  % IPS
		compiled_roi(:, hemi, subj, 14) = (roi == 25);              % FEF
		% Now eliminate labels where incorrect tval dominates
		bodytval = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.flocbodiestval.mgz']);
		facetval = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.flocfacestval.mgz']);
		wordtval = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.flocwordtval.mgz']);
		placetval = cvnloadmgz([data_dir 'subj0' num2str(subj) '/fsaverage/' hemis{hemi} '.flocplacestval.mgz']);
		% Face areas
		compiled_roi((facetval < bodytval) | (facetval < placetval) | (facetval < wordtval), hemi, subj, 5) = 0;
		compiled_roi((facetval < bodytval) | (facetval < placetval) | (facetval < wordtval), hemi, subj, 6) = 0;
		compiled_roi((facetval < bodytval) | (facetval < placetval) | (facetval < wordtval), hemi, subj, 7) = 0;
		% Body areas
		compiled_roi((bodytval < facetval) | (bodytval < placetval) | (bodytval < wordtval), hemi, subj, 8) = 0;
		compiled_roi((bodytval < facetval) | (bodytval < placetval) | (bodytval < wordtval), hemi, subj, 9) = 0;
		% Place areas
		compiled_roi((placetval < facetval) | (placetval < bodytval) | (placetval < wordtval), hemi, subj, 10) = 0;
		% Word areas
		compiled_roi((wordtval < facetval) | (wordtval < bodytval) | (wordtval < placetval), hemi, subj, 11) = 0;
	end
end
save([data_dir 'group/fsaverage/compiled_cortical_rois.mat'], 'compiled_roi');
