function Net = SFStructureN(N,targetMeanDeg,m0,p0,seed)
if nargin < 2 || isempty(targetMeanDeg)
    targetMeanDeg=10;
end
if nargin < 3 || isempty(m0)
    m0=30;
end
if nargin < 4 || isempty(p0)
    p0=0.10;
end
if nargin < 5 || isempty(seed)
    seed=3;
end
rng(seed,'twister')
A=sparse(N,N);
for i=1:m0-1
    for j=i+1:m0
        if rand<p0
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end
M0=nnz(triu(A(1:m0,1:m0),1));
targetEdges=round(targetMeanDeg*N/2);
nGrow=N-m0;
growEdges=targetEdges-M0;
baseEdges=floor(growEdges/nGrow);
remEdges=growEdges-baseEdges*nGrow;
edgePlan=baseEdges*ones(1,nGrow);
if remEdges>0
    idx=randperm(nGrow,remEdges);
    edgePlan(idx)=edgePlan(idx)+1;
end
for k=m0+1:N
    nOld=k-1;
    eNow=edgePlan(k-m0);
    degOld=sum(A(1:nOld,1:nOld)~=0,2)';
    prob=(degOld+1)/sum(degOld+1);
    oldSet=weightedSample(prob,eNow);
    A(k,oldSet)=1;
    A(oldSet,k)=1;
end
A(1:N+1:end)=0;
Net.N=N;
Net.Adj=sparse(A~=0);
end

function selected = weightedSample(prob,nPick)
prob=prob(:)';
selected=zeros(1,nPick);
picked=false(1,numel(prob));
for i=1:nPick
    p=prob;
    p(picked)=0;
    p=p/sum(p);
    j=find(cumsum(p)>=rand,1,'first');
    selected(i)=j;
    picked(j)=true;
end
end
