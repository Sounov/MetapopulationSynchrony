function dydt = ODE_sys(~,y,p)
dydt = [p.R*y(1)*(1-p.c*y(1)/p.R)*(y(1)-p.mu)/(p.nu+y(1))-p.alpha*y(1)*y(2)/(p.beta+y(1));
        p.s*y(2)*(1-p.q*y(2)/(y(1)+p.eps))];
end
