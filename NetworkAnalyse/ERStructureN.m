function Net = ERStructureN(N,targetMeanDeg,sparsity,seed)
if nargin < 2 || isempty(targetMeanDeg)
    targetMeanDeg=10;
end
if nargin < 3 || isempty(sparsity)
    sparsity=0.85;
end
if nargin < 4 || isempty(seed)
    seed=4;
end
rng(seed)
pUndir=targetMeanDeg/(N-1);
pUndir=min(max(pUndir,0),1);
pEff=1-sqrt(1-pUndir);
edgeProb=min(pEff/max(1-sparsity,eps),1);
D=rand(N,N)<edgeProb;
D(1:N+1:end)=false;
idx=find(D);
if ~isempty(idx)
    nRemove=min(round(numel(idx)*sparsity),numel(idx));
    D(idx(randperm(numel(idx),nRemove)))=false;
end
A=D|D.';
A(1:N+1:end)=false;
targetEdges=min(round(N*targetMeanDeg/2),N*(N-1)/2);
curEdges=nnz(triu(A,1));
if curEdges>targetEdges
    list=find(triu(A,1));
    remove=list(randperm(numel(list),curEdges-targetEdges));
    A(remove)=false;
    A=A&A.';
elseif curEdges<targetEdges
    missing=find(triu(true(N),1)&~A);
    nAdd=min(targetEdges-curEdges,numel(missing));
    add=missing(randperm(numel(missing),nAdd));
    A(add)=true;
    A=A|A.';
end
Net.N=N;
Net.Adj=sparse(A~=0);
end
