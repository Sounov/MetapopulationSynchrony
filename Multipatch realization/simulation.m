clc; clear; close all;
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
xi=1.1;
D=0.001;
D1=D; D2=xi*D;   % Diffusion parameters

% Simulation Parameters:================================================
Np=100;      % Number of patches

grid1 = 1; grid2=Np ; 
delh = 1;
T = 100;     % End time point
delt=0.01;
X = grid1:1:grid2;
mu1 = (delt)/(delh^2); 
J = round((grid2-grid1)/delh);
dimJ = J+1; n = dimJ; 
N = round(T/delt);

%% Climate variability ======================
rng(2)  
%===================== Indentical climatic fluctuation ====================
Rstar      = 3.7;
Rend       = 1.9;
Rmean      = (Rstar+Rend)/2;
Rm=2.3;       % critical climate value 
Tend  = T;    % end time
RR    = 1.0;  % Change expected in interval 1 year
PnBin = 0.2;  % probability of switch
[time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
climate = @(tt)interp1(time,cli1,tt);
R = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;
R1 = R(0:delt:Tend);
R1vec = ones(Np,1)*R1;  % climatic fluctuations

%% Inclusion of network topology========================
[~,~,L] = laplacian([dimJ], {'P'});  % laplacian matrix for ring network

%% Initialization =============================================:
U0 = 2 + 4*(rand(dimJ,1));  V0 = 0.001 + 0.01*(rand(dimJ,1));  % initial values 
u=zeros(n,1); v=zeros(n,1); F=zeros(n,1); G=zeros(n,1); y1=zeros(n,1); y2=zeros(n,1);
U0=U0'; V0=V0';  u=U0(:); v=V0(:);
B=sparse(1:n,1:n,1,n,n); 
B1=B-mu1*D1*L; B2=B-mu1*D2*L; 
%% =================================================
for nt=1:N
    F = R1vec(:,nt).*u.*(1-c*u./R1vec(:,nt)).*(u-mu)./(nu+u)-alpha*u.*v./(beta+u);
    G = s*v.*(1-q*v./(u+eps));
    y1 = B1*u + delt*F; y2 = B2*v + delt*G;
    u=real(y1); v=real(y2);
    U_grid(:,mod(nt,T/delt)+1) = u;
    V_grid(:,mod(nt,T/delt)+1) = v;
    Time(nt) = nt*delt;
    %% 2D pattern ===================================
    REM=rem(nt*delt,T);
    if REM==0
        data = U_grid;
        U_grid(U_grid<mu) = 0;
        tt = (nt+1)*delt-T:delt:nt*delt;
    end
end

%% Plot=================
% spatiotemporal plot=====================
figure(1)
pcolor(X,tt,U_grid'); shading flat; colormap([1,1,1;jet]); clim([0 inf]); colorbar;
xlabel('Node Index','interpreter','latex','FontSize',18); ylabel('$t$','interpreter','latex','FontSize',18,'Rotation',0);
% time series plot ================
% figure(2)
% plot(tt,U_grid')

%% Climate variablity plotter:=========================
% figure(3);
% plot(0:delt:Tend, R1);

%% Video creation: ============================================
load("Climate_timeseries.mat");
R1_amp = unique(R1,"stable");
basin_R1 = cell(1,length(R1_amp));
idx_r1 = cell(1,length(R1_amp));
x1=0; x2=14;
y1=0; y2=0.05;
n1=300;
p1=linspace(x1,x2,n1);
q1=linspace(y1,y2,n1);
load('Basin_data_zoomed.mat');
basin_R1_zoomed = basin_R1;
load('Basin_data.mat');
vid = VideoWriter('sample','MPEG-4');
vid.Quality = 100;
vid.FrameRate = 20;
open(vid)
for t=1:10:length(Time)
    idx_E = find(U_grid(:,t)<mu); 
    ln = length(idx_E);
    figure(1)
    subplot 121
    idx1 = find(R1_amp == R1(t));
    M1=basin_R1{1,idx1};
    h1=contour(p1,q1,M1,1,'--','LineWidth',2,'color',[0.85,0.33,0.10]); hold on;
    txt = ['$r_1=$' num2str(R1_amp(idx1))];
    text(10,0.045,txt,'interpreter', 'latex','FontSize',18);
   
    plot(U_grid(:,t),V_grid(:,t),'o','LineWidth',1,'MarkerSize',4); 
    title(['Time=',num2str(t*delt)],'FontSize',18);
    xlabel('$N_1$','FontSize',18,'Interpreter','latex')
    ylabel('$P_1$','FontSize',18,'Interpreter','latex')
    axis([x1 x2 y1 y2])
    hold off;
    subplot 122
    idx1 = find(R1_amp == R1(t));
    M1=basin_R1_zoomed{1,idx1};
    h1=contour(linspace(0,0.3,n1),linspace(0,0.01,n1),M1,1,'--','LineWidth',2,'color',[0.85,0.33,0.10]); hold on;

    plot(U_grid(:,t),V_grid(:,t),'o','LineWidth',1,'MarkerSize',4); 
    xline(0.03,'--');
    title(['Number of extinct patch=' num2str(ln)],'FontSize',16, Interpreter='latex');
    xlabel('$N_1$','FontSize',18,'Interpreter','latex')
    ylabel('$P_1$','FontSize',18,'Interpreter','latex')
    axis([0 0.3 0 0.01])
    hold off;
    %=======================================
    left   =  .2e3;
    bottom =  0.1e3;
    width  =  1200;
    height =  .5e3;   
    set(gcf,'Position',[left bottom width height],'renderer', 'painters')
    set(gcf,'color','w');
    curframe = getframe(gcf);
    writeVideo(vid, curframe)
    hold off;
end
writeVideo(vid, curframe)
close(vid)



