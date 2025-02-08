% Parameters:
clear all;
global r c alpha beta s q mu nu eps
r=2;
c=0.22;
alpha=505;
beta=0.3;
s=0.85;
q=205;
mu=0.03;
nu=0.003;
eps=0.031;
%% Climate variability construction: ======================
Rstar      = 3.7;
Rend       = 1.9;
Rm         = 2.3;
Rmean      = (Rstar+Rend)/2;
Tend  = 100;    % end time
RR    = 1.0;      % %avreage length of Type-L/H period
PnBin = 0.2;    % probability of switch
[time,swich,cli1,cli] = Climate(Tend,RR,PnBin);
climate = @(tt)interp1(time,cli1,tt);
R = @(tt)Rmean + climate(tt)*(Rstar-Rend)/2;

% Ode simulation with climate variability: ===============================
tspan = [0 Tend];
yint = [4, 0.001];
opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'NonNegative',[1,2]);
[tsol,ysol] = ode45(@(t,y)odefun(t,y,R), tspan, yint, opts);

%% Climate variablity plotter:=========================
figure(1);
plot([0 Tend],[Rstar,Rstar],'k:','LineWidth',1); hold on
plot([0 Tend],[Rend,Rend],'k:','LineWidth',1); hold on;
plot(time,Rmean + cli1.*((Rstar-Rend)/2),':k','LineWidth',1); hold on;
for ind = 1:length(swich)
    X = linspace(sum(swich(1:ind-1)),sum(swich(1:ind)),swich(ind));
    Y = Rmean + cli(sum(swich(1:ind-1))+1:sum(swich(1:ind)))*(Rstar-Rend)/2;
    if Y(1)>=Rm
        figure(1)
        plot(X,Y,'r-','LineWidth',3)
    else
        figure(1)
        plot(X,Y,'b-','LineWidth',3)
    end
    if swich(ind) ==1
        X1 = [X-1,X];
        Y1 = [Y, Y];
        if Y1(1)>=Rm
            figure(1)
            plot(X1,Y1,'r-','LineWidth',2)
        else
            figure(1)
            plot(X1,Y1,'b-','LineWidth',2)
        end
    end
end
box on
axis([0 Tend Rend-0.3 Rstar+0.3])
%===========================
% rescaled into rmin to rmax
% figure(1)
% plot(tsol,rescale(ysol(:,1),1.9,3.7,'InputMin',0.0610,'InputMax', max(ysol(:,1))),'-','LineWidth',2,'Color',[0 0 1]); hold on;
% plot(tsol,rescale(ysol(:,2),1.9,3.7,'InputMin',5.3035e-04,'InputMax',0.0395),'-','LineWidth',2,'Color',[1 0 0]); hold on;

figure(2)
plot(tsol,ysol(:,1),'-','LineWidth',2,'Color',[0 0 1]); hold on;
plot(tsol,10^2*ysol(:,2),'-','LineWidth',2,'Color',[1 0 0]); hold on;



%% Tipping time finder ============================================
[m,n] = size(ysol);
for ii=1:1:n
    X = ~isnan(ysol(:,ii)');
    Y = cumsum(X-diff([1,X])/2);
    ysol(:,ii) = interp1(1:nnz(X),ysol(X,ii)',Y)';
end
cond1 = ysol(end,1);
if cond1<mu            % Event other than this is persistence
    tip1 = TipTime(tsol,ysol(:,1));
    disp(['Tipping time=', num2str(tip1)])
end

