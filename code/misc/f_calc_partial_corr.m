function corrmat = f_calc_partial_corr(X, y, Z)
% Calculates partial correlation between X and y given Z for pulvinar project.
%
% Note that despite labeling the matrices X and y, we will treat X and the to-be-predicted data
% and y as the predictor of interest.
%
% See: https://en.wikipedia.org/wiki/Partial_correlation
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

    % Loop through all voxels, estimating for each the partial correlation between target region y and nuisance regions (columns of Z)
    for ii = 1:n_voxel
        % Calculate coefficients for linear regression between Z and y and Z and x
        b_star_y = Z \ y;
        b_star_x = Z \ X(:, ii);

        % From these coefficients, get residuals
        e_x = X(:, ii) - Z * b_star_x;
        e_y = y - Z * b_star_y;

        % Now, find correlation coefficient between the two and square to get R2
        temp = corrcoef(e_x, e_y);
        corrmat(ii) = temp(1, 2)^2;
    end
end
