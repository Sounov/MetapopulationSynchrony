clc; clear all; close all;
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

D = [0.001; 0.002; 0.003; 0.004; 0.005; 0.006; 0.007; 0.008; 0.009; 0.01; 0.02; 0.03; 0.04; 0.05; 0.1; 0.2; 0.3; 0.4; 0.5; 1];
SIM = 50; ITER=20;
N_pers=zeros(length(D),SIM);
N_ce=zeros(length(D),SIM);
N_pe=zeros(length(D),SIM);
N_rev=zeros(length(D),SIM);


% Simulation Parameters:================================================
Np=100;      % Number of patches
grid1 = 1; grid2=Np ; delh = 1;
T = 100; Tend=T; delt=0.001;
X = grid1:1:grid2;
mu1 = (delt)/(delh^2); J = round((grid2-grid1)/delh);
dimJ = J+1; n = dimJ; N = round(T/delt);
Rstar      = 3.7;
Rend       = 1.9;
Rmean      = (Rstar+Rend)/2;
Rm=2.3;  % critical climate value
Tend  = T;    % end time
RR    = 1.0;      % %avreage length of Type-L/H period
PnBin = 0.2;    % probability of switch


%% Initialization =============================================:
for dis = 1:1:length(D)
    D1=D(dis); D2=1.1*D(dis);   % Diffusion parameters
    for sim = 1:1:SIM
        [~,~,L] = laplacian([dimJ], {'P'});
        %         L = scale_free(100,1,1);
        %         L = erdos_reyni2(100,0.002);
        tic
        for iter=1:1:ITER
            disp(['D=', num2str(D1), 'sim=', num2str(sim), 'Iter=', num2str(iter) ]);
            %% Climate variability ======================
            %===== Homogeneous or lagged or perturbed============================
            [time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
            climate = @(tt)interp1(time,cli1,tt);
            R = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
            R1 = R(0:delt:Tend);
            R1vec = ones(Np,1)*R1;
            % Heterogeneous:===================
            %             R1vec  = ind_climate(T,100,delt);

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
                    %                     figure(2)
                    %                     plot(tt,U_grid')
                end
            end
            z = extinction(U_grid(:,2:end));
            d = diff(z);
            d(d<0) = -1;
            d(d>0) = 1;
            Nrev = nnz(d<0);
            if nnz(z) == 0
                % persistence
                N_pers(dis,sim) = N_pers(dis,sim)+1;
            else
                % extinction
                if z(end)==100
                    % complete extinction
                    N_ce(dis,sim) = N_ce(dis,sim)+1;
                else
                    % partial extinction
                    N_pe(dis,sim) = N_pe(dis,sim)+1;
                end
            end
            if Nrev > 0
                N_rev(dis,sim) = N_rev(dis,sim)+1;
            end
        end
        toc
    end
end



