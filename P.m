function heatker = P(t,x)
    % Computes the heat kernel
    % for time parameter t
    % and space parameter x
    
    heatker = ( 2 * t ).^(-1/2) .* exp( - dist2(x) / (2 * t) ) ;
    
    heatker = heatker ./ sum(heatker(:)) ;

end