function M = randvar(dim,range)

    % generates random variance-covariance matrix of the specified dimesion dim
    % within the specified range.

    d = range * rand(dim,1); % The diagonal values
    t = triu(bsxfun(@min,d,d.').*rand(dim),1); % The upper trianglar random values
    M = diag(d)+t+t.'; % Put them together in a symmetric matrix
   
end