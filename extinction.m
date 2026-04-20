function numvec = extinction(vec)
%% Nonzero elements in a row=======================
[m,n]=size(vec);
numvec = zeros(n,1);
for i=1:1:n
    numvec(i,1) = m-nnz(vec(:,i));
end
end