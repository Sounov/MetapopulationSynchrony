function vec = adjust(y)
n = length(y);
for ii=1:1:n
    if isnan(y(ii,1))==1
        y(ii,1)=y(ii-1,1);
    end
end
vec = y;