function tent = tent(x,x0,epsilon)

    % Tent function for : 
    % - space grid x
    % - center in x_0
    % - width parameter epsilon
   
    tent = (1/epsilon) * (1 - abs(x-x0) / epsilon) ;
    
    tent(tent<0) = 0 ;
    
end