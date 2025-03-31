function [trial_idxs_subcortical, trial_idxs_cortical] = f_select_trials_method2(subj, n_sess)
% Randomly permutes combinations of trials for "Method 2" pulvino-cortical correlation analyses
%
% In this function, vectors of trial indices are generated that allow for correlation of 
% responses between cortical targets and subcortical targets using "same image" but 
% "different trial" responses. 
%
% Args:
%	subj (int): whichs subject are we selecting trials for
%	n_sess (int): how many sessions did this subject complete
%
% Returns:
%	trial_idxs_subcortical (array): array of shape (10000, 6). Each column consists represents 
%		an indexing into the dataset. When a given column is paired with the same column
%		from the trial_idxs_cortical array below, the two columns index into a set of 
%		"same-image different-trial" responses between subcortex and cortex.
%	trial_idxs_cortical(array): see above, cortical version

% Config
config_Guestetal2025_NSDPulvinar;

% Calculate number of trials this person actually completed
n_trial = 750*n_sess;

% Load experimental design information and extra this subject's info
exp_design = load([nsd_dir 'experiments/nsd/nsd_expdesign.mat']);
ordering = exp_design.subjectim(:, exp_design.masterordering);
img_idx_each_trial = ordering(subj, 1:30000);  % 73000-image-index for image presented on each trial to this subject
imgs_unique = unique(img_idx_each_trial);  % 73000-image-index for each unique image presented to this subject
trial_idxs = 1:30000;  % 30000-trial-index for each trial

% First we construct an empty matrix of indices (30000-trial-indices) for subcortex and cortex
trial_idxs_cortical = zeros(10000, 6);  % 10000 unique images x 6 possible permutations
trial_idxs_subcortical = zeros(10000, 6);  % 10000 unique images x 6 possible permutations

% For each image, we do the following:
%	1. The goal of this code is to arrive at two (10000, 6) matrices of 30000-trial-indices 
%		into the dataset. 
%	2. We loop over 10000-unique-image indices and then identify how many repeats this 
%		subject actually saw for a given unique image and branch:
%		2.a. If the subject only saw one repeat, skip and don't use this data
%		2.b. If the subject saw two repeats, generate both permutations of the 
%		     30000-trial-indices indicating on which trials the subject saw these 
%  		     repeats. Now we have a matrix of size (2, 2). Each row consists of a 
%		     pair of 30000-trial-indices encoding one subcortical response (first element)
%		     and one cortical response (second element). We randomize the order of the
%		     rows. Then we extract each column, append four NaNs (reflecting the possible
%		     permutations we're missing out on because the subject didn't see the last
%		     repeat), and assign this vector as the u-th row of the corresponding output
%		     matrix (trial_idxs_subcortical if it's the first column, etc.), where u
%		     indicates the 10000-unique-trial index that we're currently on in the loop.
%		2.c. If the subject saw all three repeats, the same procedure as for 2.b is taken 
%		     except that the intermediate matrix is (6, 2) and thus we don't need to
%		     add NaNs.
%	Note, the randomization of the rows in the interim columns ensures that there is
%	no bias in any estimation procedure introduced by arbitrary ordering of the 
%	30000-unique-trial indices corresponding to each 10000-unique-image index. 

% Loop through images presented to this subject
for img=1:10000
	% Get the 30000-trial-indices where this image was presented
	idxs = trial_idxs(img_idx_each_trial == imgs_unique(img));
	% Get rid of indices that are beyond the number of trials subject actually did
	idxs = idxs(idxs <= n_trial);
	% Branch based on the number of repeats
	if length(idxs) <= 1
		% If we have just one index, we can't perform same-image different-trial correlations, so just put NaNs in the matrix
		trial_idxs_subcortical(img, :) = NaN;
		trial_idxs_cortical(img, :) = NaN;
	elseif length(idxs) == 2
		% If we have two indices, we only have two permutations, so put those in the matrix and fill the rest of the spots with NaNs
		temp = perms(idxs);
		temp = temp(randperm(size(temp, 1)), :);
		trial_idxs_subcortical(img, :) = [temp(:, 1)' NaN NaN NaN NaN];
		trial_idxs_cortical(img, :) = [temp(:, 2)' NaN NaN NaN NaN];
	else length(idxs) == 3
		% If we have three indices, we have six permutations, so put those in the matrix
		temp = perms(idxs);
		temp = temp(:, 1:2);
		temp = temp(randperm(size(temp, 1)), :);
		trial_idxs_subcortical(img, :) = temp(:, 1);
		trial_idxs_cortical(img, :) = temp(:, 2);	
	end
end
