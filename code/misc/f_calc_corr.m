function corrmat = f_calc_corr(X, Y)
% Calculates correlation between two matrices of data with low-memory for loop
%
% Accepts two matrices with m and k rows respectively. Each row is the same length. All m rows of the first matrix are correlated with all k rows of the second matrix to result in (m, k) correlation values. This is done using a for loop to minimize the amount of memory that is used in the correlation calculation.

% Args:
%	X (array): array of shape (m, n)
%	Y (array): array of shape (k, n)
%
% Returns:
%	corrmat (array): array of shape (m, k)

	% Calculate the number of samples (n)
	n_samp = size(X, 2);
	% Demean X and Y
	X = X-mean(X, 2);
	Y = Y-mean(Y, 2);
	% Construct an empty array of the desired output size (m, k)
	cov_xy = zeros(size(X, 1), size(Y, 1));
	% Loop through m and k
	for ii=1:size(X, 1)
		for jj=1:size(Y, 1)
			cov_xy(ii, jj) = 1/(n_samp-1) * X(ii, :) * Y(jj, :)';
		end
	end
	%cov_xy = 1/(n_samp-1) * X*Y';
	var_x = 1/(n_samp-1) * sum(X.^2, 2);
	var_y = 1/(n_samp-1) * sum(Y.^2, 2);
	corrmat = cov_xy./(sqrt(var_x * var_y'));
end
