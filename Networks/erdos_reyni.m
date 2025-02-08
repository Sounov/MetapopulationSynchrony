function [L,A] = erdos_reyni(N,p)
% Input:=========
% N: Number of nodes
% p: Probablity of an edge between any two nodes
% Output:=======
% L: Laplacian matrix 
% A: adjacency of matrix 
cond=0; 
count=0;
Adj = zeros(N,N);
Adj = diag(ones(N-1,1),1);
Adj = Adj+Adj';
while N ~= cond
    count=count+1;
    H = Adj;
    % Loop over each possible pair of nodes
    for i=1:N
        for j=i+1:N
            % Generate a random number between 0 and 1
            r=rand(1);
            % If r is less than or equal to p, add an edge
            if r<=p
                H(i,j)=1;
                H(j,i)=1; % Because the graph is undirected
            end
        end
    end
    A=H;
    % Visualize the graph
    B = graph(A);
    bins = conncomp(B);
    cond = sum(bins);
end
L=laplacian(B);
Deg = sum(A);
figure(1)
plot(B);
p = plot(B,'Layout','force','NodeLabel',{});
p.MarkerSize=5;










