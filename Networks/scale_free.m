function [L,A]=scale_free(N,m,s)
% Input:=======
% N: total number of nodes
% s: number of connection of newly connected nodes
% m: number of initial connected nodes
% Output:=======
% L: Laplacian matrix 
% A: adjacency of matrix 

Adj = zeros(N,N);
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
B = graph(Adj);
L=laplacian(B);
figure(1)
plot(B);
p = plot(B,'Layout','force','NodeLabel',{});
p.MarkerSize=10;
