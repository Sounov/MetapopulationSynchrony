%% Parameters:
clear all;
warning off
global r c alpha beta s q mu nu eps D1 D2
r=2;
c=0.22;
alpha=505;
beta=0.3;
s=0.85;
q=205;
mu=0.03;
nu=0.003;
eps=0.031;
mult=10^2;
%% =============================
rng(365) % For Figure 4abcd
rng(251) % For Figure 4efgh
rng(2) % For Figure 4ijkl

D=0.001; D1=D; D2=1.1*D;
% index=[];

%% Climate variability construction: ======================
% for i=1:1:10000

Rstar      = 3.7;
Rend       = 1.9;
Rmean      = (Rstar+Rend)/2;
Tend  = 100;    % end time
RR    = 1;      % %avreage length of Type-L/H period
PnBin = 0.2;    % probability of switch
[time1,swich1,cli11,cli12] = Climate(Tend,RR,PnBin);
[time2,swich2,cli21,cli22] = Climate(Tend,RR,PnBin);
climate1 = @(tt)interp1(time1,cli11,tt);
climate2 = @(tt)interp1(time2,cli21,tt);

R1 = @(tt)Rmean + climate1(tt)*(Rstar-Rend)/2;
R2 = @(tt)Rmean + climate2(tt)*(Rstar-Rend)/2;
%% Ode simulation with climate variability: ===============================
tspan = [0 Tend];
yint = [4*rand(1)+3, 0.001, 4*rand(1)+3, 0.001];
opts = odeset('RelTol',1e-6,'AbsTol',1e-5,'NonNegative',[1:4]);
[tsol,ysol] = ode15s(@(t,y)odefun_two_patch(t,y,R1,R2), tspan, yint, opts);
ysol(:,1)=adjust(ysol(:,1));
ysol(:,3)=adjust(ysol(:,3));

% if (ysol(end,1)<mu) && (ysol(end,3)<mu)
% disp(['rng=', num2str(i)]);
% index = [index; i];

%% Climate variablity plotter:=========================
figure(1)
hold on
for i=1:2
    if i==1
        swich = swich1; time=time1; cli1 = cli11; cli=cli12;
    else
        swich = swich2; time=time2; cli1 = cli21; cli=cli22;
    end
    figure(1)
    subplot(2,1,i)
    plot(time,Rmean + cli1.*((Rstar-Rend)/2),':k','LineWidth',1); hold on;
    for ind = 1:length(swich)
        X = linspace(sum(swich(1:ind-1)),sum(swich(1:ind)),swich(ind));
        Y = Rmean + cli(sum(swich(1:ind-1))+1:sum(swich(1:ind)))*(Rstar-Rend)/2;
        if Y(1)>=2.3
            %             figure(i)
            plot(X,Y,'r-','LineWidth',2); hold on;
        else
            %             figure(i)
            plot(X,Y,'b-','LineWidth',2); hold on;
        end
        if swich(ind) ==1
            X1 = [X-1,X];
            Y1 = [Y, Y];
            if Y1(1)>=2.3
                %                 figure(i)
                plot(X1,Y1,'r-','LineWidth',2); hold on;
            else
                %                 figure(i)
                plot(X1,Y1,'b-','LineWidth',2); hold on;
            end
        end
    end
    axis([0 Tend Rend Rstar])
    Pmax = max(ysol(:,2*i)); Pmin = min(ysol(:,2*i));
    Nmax = max(ysol(:,2*i-1)); Nmin = min(ysol(:,2*i-1));
    Ap = (Rstar-Rend)/(Pmax-Pmin); Bp = (Rstar*Pmin-Rend*Pmax)/(Pmax-Pmin);
    An = (Rstar-Rend)/(Nmax-Nmin); Bn = (Rstar*Nmin-Rend*Nmax)/(Nmax-Nmin);
%     figure(1)
%     subplot(2,1,i)
%     plot(tsol,Ap*ysol(:,2*i)-Bp,'-','LineWidth',2,'Color',[.65 .65 .65]); hold on;
%     plot(tsol,An*ysol(:,2*i-1)-Bn,'-','LineWidth',2,'Color',[0 0 0]); hold on;
end
[m,n] = size(ysol);
for ii=1:1:n
    X = ~isnan(ysol(:,ii)');
    Y = cumsum(X-diff([1,X])/2);
    ysol(:,ii) = interp1(1:nnz(X),ysol(X,ii)',Y)';
end
figure(3)
% subplot 211
plot(tsol,ysol(:,1),'-b',tsol,ysol(:,3),'-r','LineWidth',2); hold on;
% data = [mult*ysol(:,2),mult*ysol(:,4)];
% [zR, zR_t, Phid, pdt] = synch(data);
% figure(3)
% subplot 212
% plot(tsol, zR_t);
% zR
% Phid/100
% end
% end

% dt = 0.01;
% t = 0:dt:Tend;
% C1 = R1(t'); C1 = adjust(C1);
% C2 = R2(t'); C2 = adjust(C2);
% Mn = (C1+C2)/2;
% figure(4)
% plot(t,Mn)
%  axis([0 Tend Rend Rstar])
