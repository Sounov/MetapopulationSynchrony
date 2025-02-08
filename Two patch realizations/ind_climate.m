function climate_vec = ind_climate(T,n,dt)
% clear all;
% T=100; n=100; dt=0.01;
%% Climate variability construction: ======================
Rstar      = 3.7;
Rend       = 1.9;
Rmean      = (Rstar+Rend)/2;
Tend  = T;    % end time
RR    = 1.0;      % %avreage length of Type-L/H period
PnBin = 0.2;    % probability of switch
R = cell(n,1);
t=0:dt:T;
for ii=1:1:n
    [time1,swich1,cli11,cli12] = Climate(Tend,RR,PnBin);
    climate1 = @(tt)interp1(time1,cli11,tt);
    R1 = @(tt)Rmean + climate1(tt)*(Rstar-Rend)/2;
    R{ii} = R1;
    R_vec(ii,:) = R1(t); 
end
climate_vec = R_vec;



