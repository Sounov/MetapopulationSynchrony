clc; clear all; close all;
% Pearson's correlation calculation for two-patch system
% Both identical and non-identical climate is observed (comment any one of them)
%% Define Parameters of Model: ============================================
c=0.22;
alpha=505;
beta=0.3;
s=0.85;
q=205;
mu=0.03;
nu=0.003;
eps=0.031;
xi=1.0;
% Simulation Parameters:================================================ 
Np=100;      % Number of patches
grid1 = 1; grid2=Np ; delh = 1;
T = 100; Tend=T; delt=0.01;
X = grid1:1:grid2;
mu1 = (delt)/(delh^2); J = round((grid2-grid1)/delh);
dimJ = J+1; n = dimJ; N = round(T/delt);
[~,~,L] = laplacian([dimJ], {'P'});
Rstar      = 3.7;
Rend       = 1.9;
Rmean      = (Rstar+Rend)/2;
Tend  = T;
RR    = 1.0;
PnBin = 0.2;
iter=1;
while iter<=1000
    D=rand(1);
    D1=D; D2=xi*D;   % Diffusion parameters
    % ===================== Identical climate ============================
    [time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
    climate = @(tt)interp1(time,cli1,tt);
    R = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
    R1 = R(0:delt:Tend);
    R1vec = ones(Np,1)*R1;
    % =====================Non identical climate======================
%     R1vec = ind_climate(T,Np,delt);
    %% Initialization =============================================
    U0 = 1 + 5*(rand(dimJ,1));  V0 = 0 + 0.01*(rand(dimJ,1));  % initial values
    u=zeros(n,1); v=zeros(n,1); F=zeros(n,1); G=zeros(n,1); y1=zeros(n,1); y2=zeros(n,1);
    U0=U0'; V0=V0';  u=U0(:); v=V0(:);
    B=sparse(1:n,1:n,1,n,n); B1=B-mu1*D1*L; B2=B-mu1*D2*L; 
    %% =================================================
    for nt=1:N
        F = R1vec(:,nt).*u.*(1-c*u./R1vec(:,nt)).*(u-mu)./(nu+u)-alpha*u.*v./(beta+u);
        G = s*v.*(1-q*v./(u+eps));
        y1 = B1*u + delt*F; y2 = B2*v + delt*G;
        u=real(y1); v=real(y2);
        U_grid(:,mod(nt,T/delt)+1) = u;
        V_grid(:,mod(nt,T/delt)+1) = v;
        Time(nt) = nt*delt;
    end
    if U_grid(:,end)>mu 
        disp(['iteration=',num2str(iter)]);
        p_corr(iter) = corr(U_grid(1,:)',U_grid(2,:)');
        Dval(iter) = D;
        iter = iter+1;
    end
end
plot(Dval,p_corr,'o');
