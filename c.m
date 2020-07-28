function dist = c(x,gcpower)
    % Takes as input the 1-DIMENSIONAL space grid x (vector)
    % Outputs the matrix of ground cost
       
    nspace = length(x) ;
    per = x(nspace)-x(1)+x(2) ;
    
    distper = zeros(nspace,nspace,3) ;
    distper(:,:,1) = abs(x'-x-per).^gcpower ;
    distper(:,:,2) = abs(x'-x).^gcpower ;
    distper(:,:,3) = abs(x'-x+per).^gcpower ;
    
    dist = min(distper,[],3) ;

end