function dG = dG(beta,p)

    % Gradient of the function G (dual of the congestion cost F)
    
    dG = beta.^ ( ( p / (p-1) ) - 1 ) ;
    
end