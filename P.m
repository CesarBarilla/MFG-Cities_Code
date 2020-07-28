function heatker = P(t,x)

    % Computes the 1D heat kernel
    % -> USING ITS GAUSSIAN APPROXIMATION ON THE TORUS
    % for time parameter t
    % and space parameter x
    
    heatker = ( 2 * t ).^(-1/2) .* exp( - dist2(x) / (2 * t) ) ;
    
    % Normalization
    heatker = heatker ./ sum(heatker(:)) ;

end