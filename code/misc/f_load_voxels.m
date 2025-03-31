% load_voxels.m
%
% Loads data from the public-facing betas into an array. Has two modes of operation. In 
% box mode it loads a rectangular prism of voxels defined by two opposite vertices. In
% mask mode it loads any number of voxels defined by a n_voxel x 3 array of coordinates. 
%
% Arguments:
%	participant (int): determines which participant's data is loaded, no default.
% 	n_sess (int): number of sessions to load, defaults to 1, returns error if more
%		sessions are requested than exist in the data. 
%	mode (int): 1 = box mode, 2 = mask mode
% 	coords (matrix): 
%		- Box mode: expects a 2 x 3 matrix composed of [coord_start; coord_end].
%			A rectangular prism will be defined with coord_start and 
%			coord_end as its most distant vertices. All the voxels within 
%			this prism (inclusive) will be loaded into MATLAB. coord_start 
%			values should always be smaller than coord_end values.
%		- Mask mode: expects 3D matrix consisting of ones and zeros indicating
%			which voxels should be loaded and which should not. Be careful
%			because mask must load in the smallest box containing all of the
%			1s which can be a lot of data.
%	resolution (int): 1 = 1mm, 2 = 1.8mm 
%	prep (int): 1 = subject-native, 2 = MNI
% 	type (int): which type of betas are requested (1 = assumehrf, 2 = fithrf, 3 = GLMdenoise),
%		defaults to fithrf.
%
% Returns:
% 	data (matrix): An array containing the requested data. 
%		- Box mode: Has the shape [session, x, y, z, trial].
%		- Mask mode: Has the shape [session, n_voxel, trial].
%
% Warnings:
%	- This function is incomplete and not all the options work! AFAIK it's only used
%	in practice in mask mode with func1mm data

function [data, d1, d2, d3, idxs] = f_load_voxels(participant, n_sess, type, mode, resolution, prep, coords, offset)

% Handle arguments
if nargin < 8
	offset = 0;
end

% Handle directories
config_Guestetal2025_NSDPulvinar;

% Define directories and handle type, throw error if type is not acceptable
data_dir_head = fullfile(nsdbeta_dir, 'ppdata');
if prep == 1
	if resolution == 1
		if type == 1
			data_dir_tail = fullfile('func1mm', 'betas_assumehrf');
			n_trial = 750;
		elseif type == 2
			data_dir_tail = fullfile('func1mm', 'betas_fithrf');
			n_trial = 750;
		elseif type == 3
			data_dir_tail = fullfile('func1mm', 'betas_fithrf_GLMdenoise_RR');
			n_trial = 750;
		else
			error('type must be one of three permitted types.');
		end
	elseif resolution == 2
		if type == 1
			data_dir_tail = fullfile('func1pt8mm', 'betas_assumehrf');
			n_trial = 750;
		elseif type == 2
			data_dir_tail = fullfile('func1pt8mm', 'betas_fithrf');
			n_trial = 750;
		elseif type == 3
			data_dir_tail = fullfile('func1pt8mm', 'betas_fithrf_GLMdenoise_RR');
			n_trial = 750;
		else
			error('type must be one of four permitted types.');
		end
	end
elseif prep == 2
	error('MNI preparation is not available.');
end

% Confirm that all of the requested sessions actually exist, throw error if they don't
files = dir(fullfile(data_dir_head, ['subj' sprintf('%02d', participant)], data_dir_tail, 'betas_session*.hdf5'));
assert(length(files) >= n_sess, 'Error: There are fewer sessions than you requested!')

% ///// Handle box mode
if mode == 1 
	% Check asssumption that size of coordinates is 2x3
	assert(size(coords, 1) == 2 && size(coords, 2) == 3, ...
        'Error: too many coordinates were passed for box mode');
    % Extract start point and end point from coords 
	coord_start = coords(1, :);
	coord_end = coords(2, :);
	% Check assumption that ending coordinate is larger than starting coordinate, throw error if not
	assert(~any(~(coord_end > coord_start)), ...
        'Error: At least one ending coordinate is larger than one starting coordinate!')
	% Preallocate space for data, which should be (n_sess x x_dim x y_xim x z_dim x n_trial [trials])
	data = zeros(n_sess, ...
		coord_end(1) - coord_start(1) + 1, ...
		coord_end(2) - coord_start(2) + 1, ...
		coord_end(3) - coord_start(3) + 1, ...
		n_trial);
	% Loop through sessions to load data
	for sess=1:n_sess
		% Print what we're doing to screen
		fprintf(['Loading data from session ' num2str(sess) ' and subject ' num2str(participant) '\n']);
		% Construct filename
		file_name = fullfile(data_dir_head, ['subj' sprintf('%02d', participant)], data_dir_tail, ...
				     ['betas_session' sprintf('%02d', sess + offset) '.mat']);
		% Print filename to screen
		fprintf(['Filename: ' file_name '\n']);
		% Connect to the data file using the matfile function (but don't load everything just yet)
		temp = matfile(file_name);
		% Load the data we want, then cast it to single, divide by 300, and store to get units of % signal change
		data(sess, :, :, :, :) = single(...
			temp.betas(coord_start(1):coord_end(1), ...
				coord_start(2):coord_end(2), ...
			   	coord_start(3):coord_end(3), ...
			   	 :))/300;
	end
elseif mode == 2
	% Compute brick and indices --- i.e., find the smallest brick containing all desired voxels (indicated by 1s) 
    % in the mask and then compute the vector indices corresponding to those 1s within the brick.
	[d1, d2, d3, idxs] = computebrickandindices(coords);
	% Preallocate space for data, should be (n_sess x n_voxel x n_trial [trials])
	data = zeros(n_sess, ...
		length(idxs), ...
		n_trial);
	% Loop through sessions to load data	
	for sess=1:n_sess
		% Print what we're doing to screen
		fprintf(['Loading data from session ' num2str(sess) ' and subject ' num2str(participant) '\n']);
		% Construct filename
		file_name = fullfile(data_dir_head, ['subj' sprintf('%02d', participant)], data_dir_tail, ...
				     ['betas_session' sprintf('%02d', sess + offset) '.hdf5']);
		% Print filename to screen
		fprintf(['Filename: ' file_name '\n']);
		% Connect to the data file using the matfile function (but don't load everything just yet)
		info = hdf5info(file_name);
		temp = h5read(file_name, '/betas', [d1(1), d2(1), d3(1), 1], [d1(end)-d1(1)+1, d2(end)-d2(1)+1, d3(end)-d3(1)+1, 750]);
		% Squish and extract voxels in mask
		data_as_int = subscript(squish(temp, 3), {idxs ':'});
		% Assert that the number of unique values in our data (before converting to single) is sufficient
		for ii=1:length(idxs)
			if length(unique(data_as_int(ii, :))) < 75
				fprintf(['WARNING: Number of unique values in voxel # ' num2str(ii) ' is only ' num2str(length(unique(data_as_int(ii, :)))) '\n']);
			end
		end
		% Load the data we want, then cast it to single, divide by 300, and store to get units of % signal change
		data(sess, :, :) = single(data_as_int)/300;
	end
else
	error('The mode you requested is not a mode. Please pick 1 = box or 2 = mask');
end
