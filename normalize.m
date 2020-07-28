function normalize = normalize(p)
    
    % Function that normalizes array p
    
    normalize = p/sum(p(:));
    
end