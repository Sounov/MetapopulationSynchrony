% function L=scale_free(N,m,s)
% global N m s k
N=50; 
s = 1; % s is the number of connection
m = 5; % m is the number of initial connected nodes

Adj = zeros(N,N);
% Adj(1:m,1:m)=1-eye(m);
for i=1:m
    for j=i+1:m
        Adj(i,j)=1;
        Adj(j,i)=1;
    end
end
for i=m+1:1:N
    list = randsample(i-1,s);
    for j=1:1:s
        jj = list(j);
        Adj(jj,i) = 1;
        Adj(i,jj) = 1;
    end
end
S=Adj;
% figure(1)
% imagesc(Adj);
B = graph(Adj);
L=laplacian(B);
figure(2)
plot(B);
p = plot(B,'Layout','force','NodeLabel',{});
p.MarkerSize=10;
% Deg = sum(Adj);
% mean(Deg)
% figure(10)
% plot(sort(Deg));
% figure(70)
% histogram(Deg,200);