function G = G(beta,p,a)
    
    G = ( 1 / a ) .^ ( 1 / (p-1) ) .* ( (p-1) / p ) * beta.^ ( p / (p-1) ) ;

    
end