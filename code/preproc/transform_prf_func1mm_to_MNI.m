% transform_prf_to_MNI.m
% 
% Uses nsd_mapdata to transform pRF parametre outputs from 1mm isotropic functional space to MNI space.
%
% Note: all NSD data is LPI arrangement; however, when we transform to MNI the output is RPI arragement.
% Therefore, before saving the outputs to disk with nsd_savenifti we flip them to LPI coordinates.

% Set target directory
store_dir = '/home/surly-raid3/dguest-data/Guestetal2021_data/';
% Define which pRF data we want to map
datas = {'prf_angle', 'prf_eccentricity', 'prf_R2', 'prf_size'};
avgs = {'nearest', 'linear', 'linear', 'linear'};
% Use nsd_mapdata to map each subject's pRF results to MNI space
for data_type=1:4
	combined_map = zeros(182, 218, 182, 8);
	for subj=1:8
		combined_map(:, :, :, subj) = nsd_mapdata(...
			subj, 'func1pt0', 'MNI', ...
			['/home/surly-raid4/kendrick-data/nsd/nsddata/ppdata/subj0' num2str(subj) '/func1mm/' datas{data_type} '.nii.gz'], ...
			avgs{data_type}, 0, ...
			[store_dir 'subj0' num2str(subj) '/mni/' datas{data_type} '.nii.gz']);
	end
	% Calculate average map and save
	% <<< NOTE >>> nsd_mapdata will return volumes in RPI arrangement because we requested MNI 
	% (per Kendrick 04/22/2020), for each volume we flip before passing to nsd_savenifti to save in
	% LPI arrangement
	if data_type == 1
		% If data is pRF angle, save 3 variants:
		%	- standard (combined_map_full): angular average, wrapped to 2Pi
		%	- rotated (combined_map_rotate): angular average after 90 degree counterclockwise rotation, wrapped to -Pi to Pi. Here,
		%		positive values from 180 to 0 indicate right visual field from top to bottom, negative values from 180 to 0 indicate
		%		left visual field from top to bottom.
		%	- absolute value of rotated (combined_map_abs_rotate): same as above, except absolute value 
		combined_map_full = wrapTo2Pi(angle(mean(exp(i*(combined_map)*pi/180), 4)))*180/pi;
		combined_map_rotate = wrapToPi(angle(mean(exp(i*(combined_map+90)*pi/180), 4)))*180/pi;
		combined_map_abs_rotate = abs(angle(mean(exp(i*(combined_map+90)*pi/180), 4)))*180/pi;
		nsd_savenifti(flip(combined_map_full, 1), [1 1 1], [store_dir 'group/mni/' datas{data_type} '.nii.gz'], 1, [92, 127, 73]);
		nsd_savenifti(flip(combined_map_rotate, 1), [1 1 1], [store_dir 'group/mni/' datas{data_type} '_rotate.nii.gz'], 1, [92, 127, 73]);
		nsd_savenifti(flip(combined_map_abs_rotate, 1), [1 1 1], [store_dir 'group/mni/' datas{data_type} '_abs_rotate.nii.gz'], 1, [92, 127, 73]);
		for subj=1:8
			nsd_savenifti(flip(squeeze(wrapTo2Pi(angle(exp(i*(combined_map(:, :, :, subj))*pi/180)))*180/pi), 1), [1 1 1], [store_dir 'subj0' num2str(subj) '/mni/' datas{data_type} '.nii.gz'], 1, [92, 127, 73]);
			nsd_savenifti(flip(squeeze(wrapToPi(angle(exp(i*(combined_map(:, :, :, subj)+90)*pi/180)))*180/pi), 1), [1 1 1], [store_dir 'subj0' num2str(subj) '/mni/' datas{data_type} '_rotate.nii.gz'], 1, [92, 127, 73]);
			nsd_savenifti(flip(squeeze(abs(angle(exp(i*(combined_map(:, :, :, subj)+90)*pi/180)))*180/pi), 1), [1 1 1], [store_dir 'subj0' num2str(subj) '/mni/' datas{data_type} '_abs_rotate.nii.gz'], 1, [92, 127, 73]);
		end
	else
		% If data is not pRF angle, simply take mean and save result
		combined_map = mean(combined_map, 4);
		nsd_savenifti(flip(combined_map, 1), [1 1 1], [store_dir 'group/mni/' datas{data_type} '.nii.gz'], 1, [92, 127, 73]);
	end
end

