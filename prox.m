function b_prox = prox(kertemp,b_other,init,T,N,theta)

dt = T/N ;

% Double entropic prox of kernels kertemp1 and kertemp2

    
   Grad = @(b) dt .* dF(b+b_other) + theta * log(b ./ kertemp) ;
    
   Hess = @(b) dt * d2F(b+b_other) + theta * b.^(-1)  ;
    
   b_prox = newton(Grad,Hess,init) ;   
   
% Function end
end
 