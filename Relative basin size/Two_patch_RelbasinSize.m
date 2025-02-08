clear all; clc; close all; warning off
% Relative basin size analysis of autonoumous 2 patch system
% the environmnetal conditions chosen randomly from (r_min, r_max)
%% Parameters:
global R1 R2 c alpha beta s q mu nu eps D1 D2
c=0.22;
alpha=505;
beta=0.3;
s=0.85;
q=205;
mu=0.03;
nu=0.003;
eps=0.031;
Tend=30;
xi = 1.1;
%% Observed dispersal rate values:=============
d=[0:0.001:0.012];
b=zeros(length(d),1);
r=zeros(length(d),1);
g=zeros(length(d),1);
tspan = [0 Tend];
opts = odeset('RelTol',1e-3,'AbsTol',1e-6,'NonNegative',[1:4]);
iter=10000;
for j=1:1:length(d)
    D1=d(j); D2=xi*d(j);
    for iter=1:1:iter
        if rem(iter,100)==0
            disp(vpa([D1,iter]));
        end
        R1=(3.7-1.9)*rand(1)+1.9; 
        R2=(3.7-1.9)*rand(1)+1.9;
        yint = [14*rand(1), 0.05*rand(1), 14*rand(1), 0.05*rand(1)];
        [tsol,ysol] = ode45(@(t,y)ODE_sys(t,y), tspan, yint, opts);
        %% Time series
        if ysol(end,1)>mu && ysol(end,3)>mu  % SS state
            g(j,1)=g(j,1)+1;
        end
        if ysol(end,1)<mu && ysol(end,3)<mu  % EE state
            r(j,1)=r(j,1)+1;
        end
        if (ysol(end,1)>mu && ysol(end,3)<mu) || (ysol(end,1)<mu && ysol(end,3)>mu) % SE and ES state
            b(j,1)=b(j,1)+1;
        end
            
    end
end
states = [r,b,g];
bar(0:0.001:0.012,states/iter,'stacked')
