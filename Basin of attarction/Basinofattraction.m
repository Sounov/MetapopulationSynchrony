clear all; clc; 
%% Bifurcation diagram:
par.R=2.3;  % Value of climatic drifting parameter r
par.c=0.22;
par.alpha=505;
par.beta=0.3;
par.s=0.85;
par.q=205;
par.mu=0.03;
par.nu=0.003;
par.eps=0.031;
eat = par.mu;
%% Reference domain choice:
x1=0; x2=14;
y1=0; y2=0.05;
n1=100;    % No. of grid selection (increase grids for smoothness)
M = zeros(n1,n1);   % Basin data matrix
p1=linspace(x1,x2,n1);
q1=linspace(y1,y2,n1);

tspan=[0 40];
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
ppm = ParforProgressbar(n1^2,'progressBarUpdatePeriod', 1.5);
%% search on space:==================================
parfor n=1:n1^2
    ppm.increment();
    j = rem(n,n1);
    i = 1 + floor(n/n1);
    if j==0
        j = n1;
        i = floor(n/n1);
    end
    x0=[p1(i); q1(j)];
    [~,ysol1] = ode45(@(t,y)ODE_sys(t,y,par),tspan,x0,opts);
    if ysol1(end,1)>=eat
        M(n)=1;
    end
end
delete(ppm);
figure(1); h = contour(p1,q1,M,1); % basin boundary plot  
prs = nnz(M)*100/n1^2;              % Number of initial values converged to limit-cycle
ext = (n1^2-nnz(M))*100/n1^2;       % Number of initial values converged to extinction
