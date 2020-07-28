function dG = dG(beta,p,a)

    % Gradient of the function G (dual of the congestion cost F)
    
    dG = ( 1 / a ) .^ ( 1 / (p-1) ) .* beta.^ ( ( p / (p-1) ) - 1 ) ;
    
end