function net = BuildNetwork(A)
A = sparse(abs(A));
N = size(A,1);
A(1:N+1:end) = 0;
A = sparse(A~=0);
if ~isequal(A,A.')
    A = sparse((A+A.')~=0);
    A(1:N+1:end) = 0;
end
[i,j] = find(triu(A,1));
M = numel(i);
if M == 0
    error('BuildNetwork: network has no edges.')
end
e = (1:M).';
B = sparse([e;e],[i;j],[ones(M,1);-ones(M,1)],M,N);
net.N = N;
net.Adj = A;
net.deg = full(sum(A,2));
net.edge_i = i;
net.edge_j = j;
net.M_edge = M;
net.B = B;
end
