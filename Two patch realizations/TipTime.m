function TT = TipTime(tsol, ysol)
mu = 0.03;
for tp = 1:1:length(ysol)
    if ysol(tp,1)<mu
        TT = [tsol(tp)];
        break
    end
end