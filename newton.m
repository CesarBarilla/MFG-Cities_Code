function newton = newton(Grad,Hess,init)
% Damped Newton Method
% Inputs are :
%   Grad (function handle) : differential of the function to minimize
%   Hess (function handle) : second derivative of the fct to minimize
%   init (array) : initialization value

b = init ;
tau = 1 ;
err = 1 ;
count = 0 ;

while err > 10^(-8)
    b_prev = b ;
    b = b_prev ...
                - tau * Hess(b_prev).^(-1) ...
                        .* Grad(b_prev)  ;
    while sum(b <= 0) > 0
        tau = tau/2 ; 
         b = b_prev ...
                - tau * Hess(b_prev).^(-1) ...
                        .* Grad(b_prev)  ;
    end
    err = norm(b_prev-b) ;
    tau = 1 ;
    count = count + 1 ;
    if mod(count,10000) == 0
        disp(['Newton iterations count :',num2str(count)]) ;
        disp(['Error :', num2str(err)]) ;
    end
end

newton = b ;

%disp(['Newton method converged after ',num2str(count),...
%    ' iterations (Error = ', num2str(err),')']) ;

end