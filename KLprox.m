function [Argmin,Min] = KLprox(mu,F,init)

func = @(p) KLdiv(p,mu) + F(p) ;

options = optimoptions(@fminunc,...
                'Display','off',...
                'Algorithm','quasi-newton');
            
[Argmin,Min] = fminunc(func,init,options) ;

end