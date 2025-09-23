function corrmat = f_calc_partial_corr(X, y, Z)
% Calculates *semi*-partial correlation between X and y given Z for pulvinar project.
%
% Note that despite labeling the matrices X and y, we will treat X and the to-be-predicted data
% and y as the predictor of interest.
%
% See: https://en.wikipedia.org/wiki/Partial_correlation
% See: https://personal.utdallas.edu/~herve/Abdi-PartialRegressionCoefficient2007-pretty.pdf
%
% Args:
%	X (array): array of shape (n_trial, n_voxel)
%	y (array): array of shape (n_trial, 1)
%	Z (array): array of shape (n_trial, 13)
%
% Returns:
%	corrmat (array): array of shape (n, 1)

	% Calculate the number of samples (n)
	n_samp = size(X, 1);
    n_voxel = size(X, 2);

    % Allocate storage to store output
    corrmat = zeros(n_samp, 1);
   
    % Augment Z to include an intercept (first column)
    Z = cat(2, ones(n_samp, 1), Z);

    % Predict y from Z to get the residual 
    b_star_y = Z \ y;
    e_y = y - Z * b_star_y;  % residuals of y predicted by Z
                             % we will use this vector to correlate with all columns of X

    % Loop through all voxels, estimating for each the semi-partial correlation between 
    % target region y and target data X
    parfor ii = 1:n_voxel
        % Print progress
        if mod(ii, 100) == 0
            fprintf('  Voxel %d of %d\n', ii, n_voxel);
        end

        % Find correlation coefficient between the two and square to get R2
        temp = corrcoef(X(:, ii), e_y);
        corrmat(ii) = temp(1, 2)^2;
    end
end
