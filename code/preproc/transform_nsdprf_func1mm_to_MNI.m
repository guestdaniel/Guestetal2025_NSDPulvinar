% transform_nsdprf_func1mm_to_MNI.m
%
% Maps data from Kendrick's nsdfaceprf analysis from func1mm space to MNI space and saves results to disk

% Set up directories
ppdata_dir = '/home/surly-raid4/kendrick-data/nsd/nsddata/ppdata/';
kkdata_dir = '~kendrick/ext/figurefiles/nsd/';
nsdfaceprf_dir = '/home/stone-ext4/generic/Dropbox/nsdfaceprf/freesurfer/';
output_dir = '/home/surly-raid3/dguest-data/Guestetal2021_data/'
% Configure which maps we want to transform
models = {'backgroundauto', 'bodyauto', 'contrastNEW', 'faceauto', 'foregroundauto', 'salience', 'wordauto'};
datatypes = {'angle', 'eccentricity', 'size', 'R2', 'predr'};

% Loop through and do conversions for each subject
for subj=1:8
	% Load in this subject's .mat file (which contains incides d1, d2, d3
	% indicating where data in brick belongs in the subject-native volume)
	load([kkdata_dir 'datab3nativesurface_subj0' num2str(subj) '.mat']);
	% Load in this subject's func1mm T1 for reference (LPI)
	T1 = load_untouch_nii([ppdata_dir 'subj0' num2str(subj) ...
			       '/func1mm/T1_to_func1mm.nii.gz']);
	T1 = T1.img;
	% Loop through models and datatypes
	for model=1:length(models)
	for datatype=1:length(datatypes)
		% Load this volume (LPI)
		data = load_untouch_nii([nsdfaceprf_dir 'subj0' num2str(subj) ...
					 '/label/' models{model} '_' datatypes{datatype} '.nii.gz']);
		data = data.img;
		% Subset data to only include first volum
		data = squeeze(data(:, :, :, 1));
		% Extract size of data
		data_size = size(data);
		% Create empty volume 
		vol = zeros(size(T1));
		% Insert data into empty volume
		vol(d1, d2, d3) = data;
		% Map to MNI (note, new_vol below will be RPI arrangement)
		% If 'angle' is in name, do nearest-neighbor interpolation
		if contains(datatypes{datatype}, 'angle')
			new_vol = nsd_mapdata(subj, 'func1pt0', 'MNI', vol, 'nearest');
		else 
			new_vol = nsd_mapdata(subj, 'func1pt0', 'MNI', vol, 'linear');
		end
		% Save to disk (flip beforehand to return to LPI)
		nsd_savenifti(flip(new_vol, 1), [1 1 1], ...
			      [output_dir 'subj0' num2str(subj) '/mni/' ...
			       models{model} '_' datatypes{datatype} '.nii.gz'], ...
			      1, [92, 127, 73]);	
	end
	end
end
