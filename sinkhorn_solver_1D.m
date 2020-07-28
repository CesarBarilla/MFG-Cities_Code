% Solver for the general problem 


% INPUTS (Should later become variables) 
% -> Inputs should all be given as vectorized as if model is 1D (2D handled
% outside of the solver by reshaping

clear

% Space grid 
nspace = 100 ;
xmin = 0 ;
xmax = 10 ;

% Time horizon
T = 4 ;
% Number of time steps
N = 20 ;

sigma = .5 ; % OT Regularization parameter
nu1 = .4 ; % Diffusion parameter for first population (inhabitants)
nu2 = .2 ; % Diffusion parameter for second population (firms)

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
x = linspace(xmin,xmax,nspace) ;

% Initial densities
gaussian1D = @(m,sigma) exp( -(x-m).^2/(2*sigma^2) ) ;
normalize = @(p)p/sum(p(:));

m2_0 = gaussian1D(2,.2) + gaussian1D(6,.1) + gaussian1D(9,.3) ;
m1_0 = gaussian1D(4,.1) ;

x0_1 = 3 ;
epsilon1 = 1.5 ;
x0_2 = 7 ;
epsilon2 = 1 ;

%m2_0 = (1/epsilon1) * (1 - abs(x-x0_1) / epsilon1) ;
%m1_0 = (1/epsilon2) * (1 - abs(x-x0_2) / epsilon2) ;

m1_0(m1_0 < 0) = 0 ;
m2_0(m2_0 < 0) = 0 ;

massmin = 0.001 ; 
m1_0 = normalize(m1_0+max(m1_0)*massmin) ;
m2_0 = normalize(m2_0+max(m2_0)*massmin) ;

% Time step
dt = T/N ; 

% OT Kernel
xi = exp(-c(x,gcpower)/sigma) ;

% Heat Kernels
P1 = P(nu1*dt,x) ;
%P1 = normalize(P1+max(P1)*massmin) ;

P2 = P(nu2*dt,x) ;
%P2 = normalize(P2+max(P2)*massmin) ;



% INITIALIZATION

Q1 = repmat(m1_0,[N+1,1]) ;
Q2 = repmat(m2_0,[N+1,1]) ;
 
% Initialize potentials
A1 = ones(N+1,nspace) ; 
A2 = ones(N+1,nspace) ;
V1 = ones(N+1,nspace) ;
V2 = ones(N+1,nspace) ;
    
    
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
freq_display = 50 ;
freq_plot = 500 ; % Frequence of plots during iterations
thrs = 10^(-7) ; % Error threshold

beta_sol = ones(1,nspace) ; % Initial guess for solution of minimization problem

err_Q1 = [] ;
err_Q2 = [] ;

err_Q1_temp = 1 ;
err_Q2_temp = 1 ;

count = 0 ; % Iteration count 

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
            
            kertemp1 = kertemp(k, A1.^(-dt*sigma / theta1) .* exp(V1) , P1) ;
            kertemp2 = kertemp(k, A2.^(-dt*sigma / theta2) .* exp(V2) , P2) ;
            
            funcwithgrad = @(beta) optimiterfunc(beta,p,theta1,theta2,kertemp1,kertemp2,dt) ;

            options = optimoptions(@fminunc,...
                  'Display','off',...
                  'SpecifyObjectiveGradient',true);
              
            init = beta_sol ;
 
            [beta_sol,~] = fminunc(funcwithgrad,init,options) ;
            
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
       
        subplot(nbsteps,1,1)
        plot(x,Q1(1,:),x,Q2(1,:),'linewidth',1.5) ;
        title('Evolution of densities')
        legend('Density of inhabitants','Density of firms')
        ylabel('k=0')
        
        for iplot = 2:nbsteps-1
        k_eval = floor((iplot-1)*((N+1)/plotsteps)) + 1 ;
        subplot(nbsteps,1,iplot)    
        plot(x,Q1(k_eval,:),x,Q2(k_eval,:),'linewidth',1.5) ; axis tight ;
        timetext = ['k =', num2str(k_eval-1)] ;
        ylabel(timetext)
        end 
        
        subplot(nbsteps,1,nbsteps)    
        plot(x,Q1(N,:),x,Q2(N,:),'linewidth',1.5) ; axis tight ;
        timetext = ['k = N =', num2str(N)] ;
        ylabel(timetext)
        
        drawnow
        
    end    
    
    if count > nbitermax
        break
    end
        
% Main loop end
end


%%

addpath('/Users/cesarbarilla/Documents/MATLAB/gif')
addpath('/Users/cesarbarilla/Documents/MATLAB/suplabel')

        figure
       
        subplot(nbsteps,1,1)
        plot(x,Q1(1,:),x,Q2(1,:),'linewidth',1.5) ;
        title('Evolution of densities')
        legend('Density of inhabitants','Density of firms')
        ylabel('k=0')
        
        for iplot = 2:nbsteps-1
        k_eval = floor((iplot-1)*((N+1)/plotsteps)) + 1 ;
        subplot(nbsteps,1,iplot)    
        plot(x,Q1(k_eval,:),x,Q2(k_eval,:),'linewidth',1.5) ; axis tight ;
        timetext = ['k =', num2str(k_eval-1)] ;
        ylabel(timetext)
        end 
        
        subplot(nbsteps,1,nbsteps)    
        plot(x,Q1(N,:),x,Q2(N,:),'linewidth',1.5) ; axis tight ;
        timetext = ['k = N =', num2str(N)] ;
        ylabel(timetext)
        [ax,h1]=suplabel(modelsumup);
    