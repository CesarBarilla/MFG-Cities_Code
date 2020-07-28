function prox_alt = prox_alt(init,kertemp1,kertemp2,theta1,theta2,dt,nspace)

func = @(b) dt * sum( F(b(1:nspace)+b(nspace+1:2*nspace) ) ) ...
        + theta1 * KLdiv(b(1:nspace),kertemp1) ...
        + theta2 * KLdiv(b(nspace+1:2*nspace),kertemp2) ;

options = optimoptions(@fminunc,...
                'Display','off',...
                'Algorithm','quasi-newton');

massmin = 0.0001 ;

init(init<10^-8) = init(init<10^-8) + massmin ;
            
[prox_alt,~] = fminunc(func,init,options) ;


end