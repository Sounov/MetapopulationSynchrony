function dydt = ODE_sys(~,y)
global R1 R2 c alpha beta s q mu nu eps D1 D2
dydt = [R1*y(1)*(1-c*y(1)/R1)*(y(1)-mu)/(nu+y(1))-alpha*y(1)*y(2)/(beta+y(1))+D1*(y(3)-y(1));
        s*y(2)*(1-q*y(2)/(y(1)+eps))+D2*(y(4)-y(2));
        R2*y(3)*(1-c*y(3)/R2)*(y(3)-mu)/(nu+y(3))-alpha*y(3)*y(4)/(beta+y(3))+D1*(y(1)-y(3));
        s*y(4)*(1-q*y(4)/(y(3)+eps))+D2*(y(2)-y(4))];
end
