function [func,dfunc] = optimiterfunc(beta,p,theta1,theta2,alpha1,alpha2,dt)

    % Returns the function to minimize in each iteration of the algorithm (func)
    % and its gradient (dfunc) for variable beta.
    % theta1, theta2 are the moving cost parameters in the problem
    % alpha1, alpha2 are the iteration dependent factors (projection of kernels)
    % p is the power for the congestion
    
    
    func = theta1 * alpha1 * exp( - (dt/theta1) * beta )' ...
            + theta2 * alpha2 * exp( - (dt/theta2) * beta )' ...
            + sum( G(beta,p) ) ;
                 
    dfunc = - alpha1 .* exp( - (dt/theta1) * beta ) ...
            - alpha2 .* exp( - (dt/theta2) * beta ) ...
            + dG(beta,p)  ;  
        
end