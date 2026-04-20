function [time,swich,climate,climate1] = Climate(Tend,RR,PnBin)
% Climates generates a sequence of negative binomial random variables that
% add up to Tend, has propablility of success in single trial of PnBin and
% corresponding number of successes RR
swichtemp = nbininv(rand(1,2*Tend),RR,PnBin);
indend = 1;
while (sum(swichtemp(1:indend))<Tend)
    indend = indend+1;
end
swichtemp1 = swichtemp(1:indend);
ind_sw = 1;
for ind_sw1=1:length(swichtemp1)
    if swichtemp1(ind_sw1)~=0
        swich(ind_sw) = swichtemp1(ind_sw1);
        ind_sw = ind_sw+1;
    end
end
climate1 = NaN(1,sum(swich));
ind_cl = 1;
for ind_sw = 1:length(swich)
    climatampl = rand;
    if rand<=0.5 %mod(ind_sw,2)==1 % good years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = climatampl;
            ind_cl = ind_cl + 1;
        end
    else  % bad years
        for in_par = 1:swich(ind_sw)
            climate1(ind_cl) = -climatampl;
            ind_cl = ind_cl + 1;
        end
    end
end
time = 0:0.001:sum(swich);
climate = NaN(size(time));
for ind_cl=1:length(time)-1
    climate(ind_cl) = climate1(floor(time(ind_cl)+1));
end
end