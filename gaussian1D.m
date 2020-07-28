function gaussian1D = gaussian1D(x,m,sigma) 
    
    % Gaussian function for mean m and variance sigma^2
    % Generates a Vector for the gaussian density on space grid x
    
    gaussian1D = exp( -(x-m).^2/(2*sigma^2) ) ;
    
end