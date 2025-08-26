function corrmat = f_calc_glm(X, Z)
% Predicts each of many columns of X as linear combination of few columns of Z
%
% Args:
%	X (array): array of shape (n_trial, n_voxel)
%	Z (array): array of shape (n_trial, 14)
%
% Returns:
%	corrmat (array): array of shape (n_voxel, 14)

	% Calculate the number of samples (n)
	n_samp = size(X, 1);
    n_voxel = size(X, 2);
    n_feature = size(Z, 2);

    % Allocate storage to store output
    corrmat = zeros(n_samp, n_feature+1);  % plus one for intercept
   
    % Augment Z to include an intercept (first column)
    Z = cat(2, ones(n_samp, 1), Z);

    % Loop through all voxels, estimating for each the GLM that predicts subcortical voxel
    % data as a linear combination of cortical responses X
    parfor ii = 1:n_voxel
        % Print progress
        if mod(ii, 100) == 0
            fprintf('  Voxel %d of %d\n', ii, n_voxel);
        end

        % Find correlation coefficient between the two and square to get R2
        corrmat(ii, :) = Z \ X(:, ii);
    end
end
