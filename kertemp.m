function kertemp = kertemp(k,A,ker)

% Computes integration of a potential matrix A by succesive convolution against
% the heat kernel in all dimensions except the k-th

[Nplus,~] = size(A) ;
N = Nplus - 1 ;

if k == 1
    kertemp = sum( ker .* A(N+1,:),2)' ;
    for i = 1:N-1
        kertemp = sum(ker .* A(N+1-i,:) .* kertemp,2)' ; 
    end 
end

if (k >= 2) && (k <= N)
    
    ker_backward = sum( ker .* A(N+1,:),2)' ;
    for i = 1:N-k
        ker_backward= sum(ker .* A(N+1-i,:) .* ker_backward,2)' ; 
    end 
    
    ker_forward = sum( ker' .* A(1,:),2)' ;
    for i = 2:k-1
        ker_forward = sum(ker' .* A(i,:) .* ker_forward,2)' ; 
    end 
    
    kertemp = ker_forward.*ker_backward ;
end

if k == N+1
    kertemp = sum( ker' .* A(1,:),2)' ;
    for i = 2:N
        kertemp = sum(ker' .* A(i,:) .* kertemp,2)' ; 
    end 
end

end