function KLdiv = KLdiv(p,q)

if size(p) == size(q)
    
    int = p .* ( log(p ./ q) - 1 ) + q ;
    KLdiv = sum(int(:)) ;
    
else
    
    msg = 'Dimensions of p and q must match' ;
    error(msg)
    
end

end