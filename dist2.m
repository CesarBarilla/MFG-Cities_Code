function dist = dist2(x)
    % Takes as input the 1-DIMENSIONAL space grid x (vector)
    % Outputs the matrix of QUADRATIC ground cost on the torus
    
    nspace = length(x) ;
    per = x(nspace)-x(1)+x(2) ;
    
    distper = zeros(nspace,nspace,3) ;
    distper(:,:,1) = (x'-x-per).^2 ;
    distper(:,:,2) = (x'-x).^2 ;
    distper(:,:,3) = (x'-x+per).^2 ;
    
    dist = min(distper,[],3) ;

end