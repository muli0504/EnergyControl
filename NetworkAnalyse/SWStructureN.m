function Net = SWStructureN(N,K,p,seed)
if nargin < 2 || isempty(K)
    K = 10;
end
if nargin < 3 || isempty(p)
    p = 0.05;
end
if nargin < 4 || isempty(seed)
    seed = 1;
end
assert(N>K)
assert(K>=2 && mod(K,2)==0)
assert(p>=0 && p<=1)
rng(seed)
R = K/2;
A = sparse(N,N);
edgeList = [];
for i = 1:N
    for d = 1:R
        j = mod(i-1+d,N)+1;
        if i~=j && A(i,j)==0
            A(i,j)=1;
            A(j,i)=1;
            edgeList=[edgeList;i,j];
        end
    end
end
for e = 1:size(edgeList,1)
    i = edgeList(e,1);
    j = edgeList(e,2);
    if rand<=p
        A(i,j)=0;
        A(j,i)=0;
        candidates=setdiff(1:N,[i,find(A(i,:)~=0)]);
        if ~isempty(candidates)
            new_j=candidates(randi(numel(candidates)));
            A(i,new_j)=1;
            A(new_j,i)=1;
        else
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end
A(1:N+1:end)=0;
A=sparse(A~=0);
Net.N=N;
Net.Adj=A;
end
