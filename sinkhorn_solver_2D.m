% Solver for the general problem 


% INPUTS (Should later become variables) 
% -> Inputs should all be given as vectorized as if model is 1D (2D handled
% outside of the solver by reshaping

clear

% Space grid 
nspace = 20 ;
nspace2 = nspace.^2 ;

x1min = 0 ;
x1max = 10 ;
x2min = 0 ;
x2max = 10 ;

% Time horizon
T = 2 ;
% Number of time steps
N = 8 ;

sigma = .2 ; % OT Regularization parameter
nu1 = .5 ; % Diffusion parameter for first population (inhabitants)
nu2 = .3 ; % Diffusion parameter for second population (firms)

% Moving cost parameters
theta1 = 1 ;    
theta2 = 4 ;

% Ground cost : linear, sqrt, or quadratic
groundcost = 'linear' ; 

if strcmp(groundcost,'sqrt') == 1
    gcpower = 1/2 ;
elseif strcmp(groundcost,'linear') == 1
    gcpower = 1 ;
elseif strcmp(groundcost,'quadratic') == 1
    gcpower = 2 ;
end

% Congestion power
p = 4 ;



% INTERNAL DEFINITIONS

% Space grid
x1 = linspace(x1min,x1max,nspace) ;
x2 = linspace(x2min,x2max,nspace) ; 

[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)] ;

gaussian2D = @(mu,sigma) mvnpdf(X,mu,sigma) ;
normalize = @(p)p/sum(p(:));


rng(1)
Sig1 = randvar(2,1) ; % Generates random variance-covariance matrix
Sig2 = randvar(2,1) ; % Generates random variance-covariance matrix
Sig3 = randvar(2,1) ; % Generates random variance-covariance matrix
Sig4 = randvar(2,1) ; % Generates random variance-covariance matrix

m1_0 = gaussian2D([3,1],Sig1) ...
        + gaussian2D([7,8],Sig2) ...
        + gaussian2D([2,2],Sig3) ;
m2_0 = gaussian2D([5,5],Sig4) ;

massmin = .001 ; 
m1_0 = normalize(m1_0+max(m1_0)*massmin) ;
m2_0 = normalize(m2_0+max(m2_0)*massmin) ;

m1_0 = m1_0' ;
m2_0 = m2_0' ; 

%m1_0_line = m1_0(:)' ;
%m2_0_line = m2_0(:)' ;

% Time step
dt = T/N ; 

% OT Kernel IN MATRIX FORM (size nspace^2 x nspace^2)
groundcost = c(X(:,1)',1) + c(X(:,2)',1) ;
xi = exp(-groundcost/sigma) ;

% Heat Kernels
P1 = P(nu1*dt,X(:,1)') .* P(nu1*dt,X(:,2)') ;
P1 = normalize(P1+max(P1)*massmin) ;

P2 = P(nu2*dt,X(:,1)') .* P(nu2*dt,X(:,2)') ;
P2 = normalize(P2+max(P2)*massmin) ;


% INITIALIZATION

Q1 = repmat(m1_0,[N+1,1]) ;
Q2 = repmat(m2_0,[N+1,1]) ;
 
% Initialize potentials
A1 = ones(N+1,nspace2) ; 
A2 = ones(N+1,nspace2) ;
V1 = ones(N+1,nspace2) ;
V2 = ones(N+1,nspace2) ;
    
    
nbsteps = 7 ;
plotsteps = 7 ;

parameterstext = ['T = ', num2str(T), ' ; N = ', num2str(N),...
    ' ; \theta_1 = ', num2str(theta1),...
    ' ; \theta_2 = ', num2str(theta2),...
    ' ; \sigma = ', num2str(sigma),...
    ' ; \nu_1 = ', num2str(nu1),...
    ' ; \nu_2 = ', num2str(nu2),...
    ' ; p = ', num2str(p),...
    ' ; ground cost : ', groundcost] ;

modelsumup = { parameterstext } ;
    
%% SINKHORN ITERATIONS 

nbitermax = 2000 ; % Set number of iteration
freq_display = 10 ;
freq_plot = 500 ; % Frequence of plots during iterations
thrs = 10^(-8) ; % Error threshold

beta_sol = ones(1,nspace2) ; % Initial guess for solution of minimization problem

err_Q1 = [] ;
err_Q2 = [] ;

err_Q1_temp = 1 ;
err_Q2_temp = 1 ;

count = 0 ; % Iteration count 

tic

while (err_Q1_temp > thrs) || (err_Q2_temp > thrs)
    
    count = count + 1 ;
    
    if mod(count,freq_display) == 0
        disp(['Iteration count : ', num2str(count)])
        disp(['Error (L2 norm on densities) = ', ...
            num2str(err_Q1_temp),' (Q1), ', ...
            num2str(err_Q2_temp),' (Q2)'])
    end
    
    % Store previous iterations
    Q1_prev = Q1 ;
    Q2_prev = Q2 ;
         
    % Update of a_k,b_k,c_k,d_k for k=1,...,N
    
       for k = 1:N+1
            
            A1(k,:) = ( ...
                          (  ...
                            exp(V1(k,:)) ...
                            .* ...
                            kertemp(k , A1.^(-dt*sigma / theta1) .* exp(V1) , P1) ...
                          ) ...
                         ./...
                          ( ...
                             A2(k,:) ...                            
                             * ...
                             xi'...
                           )  ...
                       )...
                       .^ ( theta1 / (theta1 + dt * sigma) ) ...
                     ;

         
            A2(k,:) = ( ...
                          (  ...
                            exp(V2(k,:)) ...
                            .* ...
                            kertemp(k , A2.^(-dt*sigma / theta2) .* exp(V2) , P2) ...
                          ) ...
                         ./...
                          ( ...
                             A1(k,:) ...                            
                             * ...
                             xi ...
                           )  ...
                       )...
                       .^ ( theta2 / (theta2 + dt * sigma) ) ...
                     ;
                 
       end
        
       V1(1,:) = ...
                log( ...
                 m1_0 ...
                 ./ ...
                 (... 
                    A1(1,:).^( - (dt*sigma) / theta1 ) ...
                     .* kertemp(1,A1.^(-dt*sigma / theta1) .* exp(V1) , P1) ...
                  ) ...
                ) ;

       V2(1,:) = ...
                log( ...
                 m2_0 ...
                 ./ ...
                  (... 
                    A2(1,:).^( - (dt*sigma) / theta2 ) ...
                     .* kertemp(1,A2.^(-dt*sigma / theta2) .* exp(V2) , P2) ...
                   ) ...
               ) ;  
       
    
        for k = 2:N+1
            
            alpha1 = kertemp(k, A1.^(-dt*sigma / theta1) .* exp(V1) , P1) ;
            alpha2 = kertemp(k, A2.^(-dt*sigma / theta2) .* exp(V2) , P2) ;
                        
            funcwithgrad = @(beta) optimiterfunc(beta,p,theta1,theta2,alpha1,alpha2,dt) ;

            options = optimoptions(@fminunc,...
                  'Display','off',...
                  'SpecifyObjectiveGradient',true);
              
            init = beta_sol ;
 
            [beta_sol,~] = fminunc(funcwithgrad,init,options) ;
                
%           init = beta_sol ;      
           
%           beta_sol = fsolve(dfunc,init,optimoptions('fsolve','Display','off')) ;
            
            V1(k,:) = (dt/theta1) * beta_sol ;
            V2(k,:) = (dt/theta2) * beta_sol ;
            
        end

    
    % Update Q1 :
        for k = 1:N+1
            Q1(k,:) = A1(k,:).^(-(dt*sigma)/theta1) ...
                         .* exp(V1(k,:)) ...
                         .* kertemp(k,A1.^(-(dt*sigma)/theta1) .* exp(V1),P1) ;
        end
    
    % Update Q here :
    
        for k = 1:N+1
            Q2(k,:) = A2(k,:).^(-(dt*sigma)/theta2) ...
                         .* exp(V2(k,:)) ...
                         .* kertemp(k,A2.^(-(dt*sigma)/theta2) .* exp(V2),P2) ;
        end
        
    % Convergence :
     err_Q1_temp = norm(Q1_prev - Q1) ;
     err_Q2_temp = norm(Q2_prev - Q2) ;
   
        
    % Plots
    
    if mod(count,freq_plot) == 0
        
    figure
        
        subplot(nbsteps,2,1)
        surf(X1,X2,reshape(Q1(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        title('Density of inhabitants') ;
        ylabel('k = 0')
        subplot(nbsteps,2,2)
        surf(X1,X2,reshape(Q2(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        title('Density of firms') ;
        
        for iplot = 2:plotsteps
        k_eval = floor((iplot-1)*((N+1)/plotsteps)) + 1 ;
        j = 2*iplot-1 ;
        subplot(nbsteps,2,j)
        surf(X1,X2,reshape(Q1(k_eval,:),[nspace,nspace]),...
            'edgecolor','none') ;
        timetext = ['k = ', num2str(k_eval-1)] ;
        ylabel(timetext)
        subplot(nbsteps,2,j+1)
        surf(X1,X2,reshape(Q1(k_eval,:),[nspace,nspace]),...
            'edgecolor','none') ;
        end
        
        subplot(nbsteps,2,2*nbsteps-1)
        surf(X1,X2,reshape(Q1(N+1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        timetext = ['k = N =', num2str(N)] ;
        ylabel(timetext)
        subplot(nbsteps,2,2*nbsteps)
        surf(X1,X2,reshape(Q2(N+1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
     
        %[ax,h1]=suplabel(modelsumup);
        
        drawnow
        
    end    
    
    if count > nbitermax
        break
    end
        
% Main loop end
end

if count < nbitermax
    disp(['Main loop converged after ',num2str(count),...
        ' iterations (Error = ', num2str(err_P_temp),...
        ' (Q1) ; ', num2str(err_Q_temp),' (Q2))']) ;
end

if count >= nbitermax
    disp(['Main loop ended after reaching max of ',num2str(count),...
        ' iterations (Error = ', num2str(err_P_temp),...
        ' (Q1) ; ', num2str(err_Q_temp),' (Q2))']) ;
end


toc

%%

addpath('/Users/cesarbarilla/Documents/MATLAB/gif')
addpath('/Users/cesarbarilla/Documents/MATLAB/suplabel')


clf

        figure
        subplot(1,2,1)
        surf(X1,X2,reshape(Q1(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ; 
        axis tight ;
        pbaspect([1,1,1])
        title('Density of inhabitants') ;
        ylabel('k = 0')
        subplot(1,2,2)
        surf(X1,X2,reshape(Q2(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ; 
        axis tight ;
        pbaspect([1,1,1])
        title('Density of firms') ;
        [ax,h1]=suplabel(modelsumup);
        gif('G2D_10.gif','Delaytime', 3/4,'frame',gcf)
        
       
        for iplot = 2:N+1
        subplot(1,2,1)
        surf(X1,X2,reshape(Q1(iplot,:),[nspace,nspace]),...
            'edgecolor', 'none') ; 
        axis tight ;
        pbaspect([1,1,1])
        title('Density of inhabitants') ;
        ylabel(['k =', num2str(iplot-1)])
        subplot(1,2,2)
        surf(X1,X2,reshape(Q2(iplot,:),[nspace,nspace]),...
        'edgecolor', 'none') ; axis tight ;
        pbaspect([1,1,1])
        title('Density of firms') ;
        [ax,h1]=suplabel(modelsumup);
        gif
        end

%%

addpath('/Users/cesarbarilla/Documents/MATLAB/suplabel')


figure
        
        subplot(nbsteps,2,1)
        contourf(X1,X2,reshape(Q1(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        title('Density of inhabitants') ;
        ylabel('k = 0')
        subplot(nbsteps,2,2)
        contourf(X1,X2,reshape(Q2(1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        colorbar
        title('Density of firms') ;
        
        for iplot = 2:plotsteps
        k_eval = floor((iplot-1)*((N+1)/plotsteps)) + 1 ;
        j = 2*iplot-1 ;
        subplot(nbsteps,2,j)
        contourf(X1,X2,reshape(Q1(k_eval,:),[nspace,nspace]),...
            'edgecolor','none') ;
        timetext = ['k = ', num2str(k_eval-1)] ;
        ylabel(timetext)
        subplot(nbsteps,2,j+1)
        contourf(X1,X2,reshape(Q2(k_eval,:),[nspace,nspace]),...
            'edgecolor','none') ;
        colorbar
        end
        
        subplot(nbsteps,2,2*nbsteps-1)
        contourf(X1,X2,reshape(Q1(N+1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        timetext = ['k = N =', num2str(N)] ;
        ylabel(timetext)
        subplot(nbsteps,2,2*nbsteps)
        contourf(X1,X2,reshape(Q2(N+1,:),[nspace,nspace]),...
            'edgecolor', 'none') ;
        colorbar
     
        [ax,h1]=suplabel(modelsumup);