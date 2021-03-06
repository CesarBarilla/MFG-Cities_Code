%% SOLVER FOR 1D Problem


% INPUTS/PARAMETERS

clear
cd('/Users/cesarbarilla/Documents/Work/Projects/MFG-Cities/MFG-Cities_Code') 

% Space grid 
    nspace = 100 ;
    xmin = 0 ;
    xmax = 10 ;
    x = linspace(xmin,xmax,nspace) ;

% Time horizon
    T = 40 ;
% Number of time steps
    N = 80 ;
% OT Regularization parameter
    sigma = .1 ; 
% Diffusion parameter for first population (inhabitants)    
    nu1 = .01 ; 
 % Diffusion parameter for second population (firms)    
    nu2 = .01 ;

% Moving cost parameters
    theta1 = 150 ;
    theta2 = 200 ;

% Ground cost : linear, sqrt, or quadratic
    groundcosttxt = 'sqrt' ; 

    if strcmp(groundcosttxt,'sqrt') == 1
        gcpower = 1/2 ;
    elseif strcmp(groundcosttxt,'linear') == 1
        gcpower = 1 ;
    elseif strcmp(groundcosttxt,'quadratic') == 1
        gcpower = 2 ;
    end

% Congestion function parameters 
    p = 50 ; % power
    a = 100 ; % multiplicative constant


% Initial conditions (comment/uncomment/modify as desired)

    %m2_0 = gaussian1D(x,2,.2) + gaussian1D(x,6,.1) + gaussian1D(x,9,.3) ;
    %m1_0 = gaussian1D(x,4,.1) ;
    
    m2_0 = gaussian1D(x,2,.2) + gaussian1D(x,5,.5) + gaussian1D(x,6,.1) + gaussian1D(x,9,.3) ;
    m1_0 = ones(1,nspace) ;

    %m1_0 = gaussian1D(x,1,.1) + gaussian1D(x,3,.1) + gaussian1D(x,5,.1) + gaussian1D(x,7,.1) ;
    %m2_0 = gaussian1D(x,2,.1) + gaussian1D(x,4,.1) + gaussian1D(x,6,.1) + gaussian1D(x,8,.1) ;

    %m1_0 = tent(x,0,2) + tent(x,10,2) ;
    %m2_0 = tent(x,5,1) ;


% Normalize initial conditions
    massmin = 0.01 ; 
    m1_0 = normalize(m1_0+max(m1_0)*massmin) ;
    m2_0 = normalize(m2_0+max(m2_0)*massmin) ;

% Time step
    dt = T/N ; 

% OT Kernel
    xi = exp(-c(x,gcpower)/sigma) ;

% Heat Kernels
    P1 = P(nu1*dt,x) ;
    P2 = P(nu2*dt,x) ;



% INITIALIZATION
    Q1 = repmat(m1_0,[N+1,1]) ;
    Q2 = repmat(m2_0,[N+1,1]) ;
    A1 = ones(N+1,nspace) ; 
    A2 = ones(N+1,nspace) ;
    V1 = ones(N+1,nspace) ;
    V2 = ones(N+1,nspace) ;


% Plotting parameters    
    nbsteps = 7 ;
    plotsteps = 7 ;

% Store parameters in text to append to graph outputs
    parameterstext = ['T = ', num2str(T), ' ; N = ', num2str(N),...
        ' ; \theta_1 = ', num2str(theta1),...
        ' ; \theta_2 = ', num2str(theta2),...
        ' ; \sigma = ', num2str(sigma),...
        ' ; \nu_1 = ', num2str(nu1),...
        ' ; \nu_2 = ', num2str(nu2),...
        ' ; p = ', num2str(p),...
        ' ; a = ', num2str(a),...
        ' ; ground cost : ', groundcosttxt] ;

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
            
            kertemp1 = kertemp(k, A1.^(-dt*sigma / theta1) .* exp(V1) , P1) ;
            kertemp2 = kertemp(k, A2.^(-dt*sigma / theta2) .* exp(V2) , P2) ;
            
            funcwithgrad = @(beta) optimiterfunc(beta,p,a,theta1,theta2,kertemp1,kertemp2,dt) ;

            options = optimoptions(@fminunc,...
                  'Display','off',...
                  'SpecifyObjectiveGradient',true);
              
            init = beta_sol ;
 
            [beta_sol,~] = fminunc(funcwithgrad,init,options) ;
            
            V1(k,:) = - (dt/theta1) * beta_sol ;
            V2(k,:) = - (dt/theta2) * beta_sol ;
            
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

if count < nbitermax
    disp(['Main loop converged after ',num2str(count),...
        ' iterations (Error = ', num2str(err_Q1_temp),...
        ' (Q1) ; ', num2str(err_Q2_temp),' (Q2))']) ;
end

if count >= nbitermax
    disp(['Main loop ended after reaching max of ',num2str(count),...
        ' iterations (Error = ', num2str(err_Q1_temp),...
        ' (Q1) ; ', num2str(err_Q2_temp),' (Q2))']) ;
end

toc

%% Output graph

addpath('/Users/cesarbarilla/Documents/MATLAB/suplabel')
cd('/Users/cesarbarilla/Documents/Work/Projects/MFG-Cities/Simulations')

        figure
       
        subplot(nbsteps,1,1)
        hold on
        plot(x,Q1(1,:),'linewidth',1.5) ;  
        plot(x,Q2(1,:),'--','linewidth',1.5) ; 
        axis tight ;
        box on ;
        hold off
        title('Evolution of densities')
        legend('Density of inhabitants','Density of firms','Location','northeast')
        ylabel('k=0')
        
        for iplot = 2:nbsteps-1
        k_eval = floor((iplot-1)*((N+1)/plotsteps)) + 1 ;
        subplot(nbsteps,1,iplot)    
        hold on
        plot(x,Q1(k_eval,:),'linewidth',1.5) ; 
        plot(x,Q2(k_eval,:),'--','linewidth',1.5) ; 
        axis tight ;
        box on ;
        hold off
        timetext = ['k =', num2str(k_eval-1)] ;
        ylabel(timetext)
        end 
        
        subplot(nbsteps,1,nbsteps)    
        hold on
        plot(x,Q1(N,:),'linewidth',1.5) ; 
        plot(x,Q2(N,:),'--','linewidth',1.5) ; 
        axis tight ;
        box on ;
        hold off
        timetext = ['k = N =', num2str(N)] ;
        ylabel(timetext)
        [ax,h1]=suplabel(modelsumup);
               
%% Output Gif

addpath('/Users/cesarbarilla/Documents/MATLAB/gif')
addpath('/Users/cesarbarilla/Documents/MATLAB/suplabel')
cd('/Users/cesarbarilla/Documents/Work/Projects/MFG-Cities/Simulations')

        figure
        
        plot(x,Q1(1,:),x,Q2(1,:),'linewidth',1) ; axis tight ;
        pbaspect([2,1,1])
        legend('Density of inhabitants','Density of firms')
        gif('Simu1D_30 .gif','Delaytime', 3/4,'frame',gcf)
        timetext = ['k =', num2str(0)] ;
        ylabel(timetext)
        xlabel(modelsumup)
        
        for iplot = 2:N+1
        plot(x,Q1(iplot,:),x,Q2(iplot,:),'linewidth',1) ; axis tight ;
        pbaspect([2,1,1])
        legend('Density of inhabitants','Density of firms')
        timetext = ['k =', num2str(iplot-1)] ;
        ylabel(timetext)
        xlabel(modelsumup)
        gif
        end
    