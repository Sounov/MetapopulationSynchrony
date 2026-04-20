clc; clear all; close all;
N_pers=0; N_ce=0; N_pe=0; N_rev=0;
%% Define Parameters of Model: ============================================
c=0.22;
alpha=505;
beta=0.3;
s=0.85;
q=205;
mu=0.03;
nu=0.003;
eps=0.031;
mult=10^2;
D12=0; D21=0;

% Simulation Parameters:================================================
D=1;
D1=D; D2=D;   % Diffusion parameters
Np=100;      % Number of patches

grid1 = 1; grid2=Np ; delh = 1;
T = 100; Tend=T; delt=0.001;
X = grid1:1:grid2;
mu1 = (delt)/(delh^2); J = round((grid2-grid1)/delh);
dimJ = J+1; n = dimJ; N = round(T/delt);

%% Initialization =============================================:
[~,~,L] = laplacian([dimJ], {'P'});

for iter=1:1:1000
    disp(['Iter=', num2str(iter)]);
    %% Climate variability ======================
    % rng(2)
    %===== Homogeneous or lagged or perturbed============================
    Rstar      = 3.7;
    Rend       = 1.9;
    Rmean      = (Rstar+Rend)/2;
    Rm=2.3;  % critical climate value
    Tend  = T;    % end time
    RR    = 1.0;      % %avreage length of Type-L/H period
    PnBin = 0.2;    % probability of switch
    [time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
    climate = @(tt)interp1(time,cli1,tt);
    R = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
    R1 = R(0:delt:Tend);
    R1vec = ones(Np,1)*R1;

    U0 = 2 + 4*(rand(dimJ,1));  V0 = 0.001 + 0.01*(rand(dimJ,1));  % initial values
    u=zeros(n,1); v=zeros(n,1); F=zeros(n,1); G=zeros(n,1); y1=zeros(n,1); y2=zeros(n,1);
    U0=U0'; V0=V0';  u=U0(:); v=V0(:);
    B=sparse(1:n,1:n,1,n,n); B1=B-mu1*D1*L; B2=B-mu1*D2*L; B12=-mu1*D12*L; B21=-mu1*D21*L;
    %% =================================================
    for nt=1:N
        F = R1vec(:,nt).*u.*(1-c*u./R1vec(:,nt)).*(u-mu)./(nu+u)-alpha*u.*v./(beta+u);
        G = s*v.*(1-q*v./(u+eps));
        y1 = B1*u + B12*v + delt*F; y2 = B21*u + B2*v + delt*G;
        u=real(y1); v=real(y2);
        U_grid(:,mod(nt,T/delt)+1) = u;
        Time(nt) = nt*delt;
        %% 2D pattern ===================================
        REM=rem(nt*delt,T);
        if REM==0
            U_grid(U_grid<mu) = 0;
            tt = (nt+1)*delt-T:delt:nt*delt;
            %             disp(['T=',num2str(nt*delt)])
            %             figure(10)
            %             pcolor(X,tt,U_grid'); shading flat; colormap([1,1,1;jet]); clim([0 inf]); colorbar;
            figure(2)
            plot(tt,U_grid')
            %             yline(mu)
            %             drawnow
        end
    end
    %     figure(10)
    %     xlabel('Node Index','interpreter','latex','FontSize',18);
    %     ylabel('$t$','interpreter','latex','FontSize',18,'Rotation',0);

    %% Climate variablity plotter:=========================
    %     figure(1);
    %     plot([0 Tend],[Rstar,Rstar],'k:','LineWidth',1); hold on
    %     plot([0 Tend],[Rend,Rend],'k:','LineWidth',1); hold on;
    %     plot(time,Rmean + cli1.*((Rstar-Rend)/2),':k','LineWidth',1); hold on;
    %     for ind = 1:length(swich)
    %         X = linspace(sum(swich(1:ind-1)),sum(swich(1:ind)),swich(ind));
    %         Y = Rmean + cli(sum(swich(1:ind-1))+1:sum(swich(1:ind)))*(Rstar-Rend)/2;
    %         if Y(1)>=Rm
    %             figure(1)
    %             plot(X,Y,'r-','LineWidth',3)
    %         else
    %             figure(1)
    %             plot(X,Y,'b-','LineWidth',3)
    %         end
    %         if swich(ind) ==1
    %             X1 = [X-1,X];
    %             Y1 = [Y, Y];
    %             if Y1(1)>=Rm
    %                 figure(1)
    %                 plot(X1,Y1,'r-','LineWidth',2)
    %             else
    %                 figure(1)
    %                 plot(X1,Y1,'b-','LineWidth',2)
    %             end
    %         end
    %     end
    %     box on
    %     axis([0 Tend Rend-0.3 Rstar+0.3])
    %=======================================================================
    z = extinction(U_grid(:,2:end));
    d = diff(z);
    d(d<0) = -1;
    d(d>0) = 1;
    Nrev = nnz(d<0);
    N_ext = nnz(d>0);
    if nnz(z) == 0
        % persistence
        N_pers = N_pers+1;
    else
        % extinction
        if z(end)==100
            % complete extinction
            N_ce = N_ce+1;
        else
            % partial extinction
            N_pe = N_pe+1;
        end
    end
    if Nrev > 0
        N_rev = N_rev+1;
    end
end
