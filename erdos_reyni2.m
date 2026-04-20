function L = erdos_reyni2(N,p)
% N=50;
% p=0.0;

% % Inputs:
% n - number of nodes
% p - probablity of an edge between any two nodes
% output:
% G- adjacency matrix of the generated graph
% Initiaize an n by n adjacency matrix with zeros
cond=0; count=0;
Adj = zeros(N,N);
Adj = diag(ones(N-1,1),1);
Adj = Adj+Adj';
Adj(end,1)=1; Adj(1,end)=1;
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
% figure(2)
% plot(B);
% p = plot(B,'Layout','force','NodeLabel',{});
% p.MarkerSize=10;
% p.NodeLabel = 0;
% p.NodeCData = rand(100,1);
% deg = sum(A);
% mean(deg)










