function Fig4
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
Gc = (0.01:0.01:0.40).';
rhoList = [0.05 0.2 0.35];
fixedGc = 0.23;
targetRho = 0.2;
recordStep = 50;
Net = SWStructureN(N,K,p,1);
net = BuildNetwork(Net.Adj);
trackIdx = unique(round(linspace(1,net.M_edge,3))).';
if numel(trackIdx) < 3
    trackIdx = (1:min(3,net.M_edge)).';
end
trackPairs = [net.edge_i(trackIdx),net.edge_j(trackIdx)];
singlePair = trackPairs(1,:);
params0 = struct('h',0.01,'T',200001,'a',0.3,'b',1.0,'c',0.1,'apha',0.1,'beta',0.77,'k1',1.0,'k2',1.8);
params0.N = N;
params0.Adj = net.Adj;
params0.deg = net.deg;
params0.edge_i = net.edge_i;
params0.edge_j = net.edge_j;
params0.M_edge = net.M_edge;
params0.B = net.B;
params0.Amp = 0.828*ones(N,1);
params0.freq = 0.19*ones(N,1);
params0.transient = round(0.3*params0.T);
params0.q0 = ones(net.M_edge,1);
params0.qMin = 1e-12;
params0.couplingType = 'energyExp';
params0.epsE = 1e-5;
params0.calcHSF = false;
params0.calcCSF = false;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = false;
params0.calcGStat = true;
params0.calcEnergySyncTime = false;
params0.calcFig7Map = false;
rng(1,'twister')
x0 = 2*rand(N,3);
if isempty(gcp('nocreate'))
    parpool;
end
numG = numel(Gc);
numRho = numel(rhoList);
CVFlat = NaN(numG*numRho,1);
parfor job = 1:numG*numRho
    [gIdx,rIdx] = ind2sub([numG,numRho],job);
    params = params0;
    params.gc = Gc(gIdx);
    params.rho_q = rhoList(rIdx);
    R = FHNres(x0,params);
    CVFlat(job) = R.gEdge_cv_time;
end
gEdgeCV = reshape(CVFlat,numG,numRho);
singleEdgeTime = cell(numRho,1);
singleEdgeG = cell(numRho,1);
for rIdx = 1:numRho
    params = params0;
    params.gc = fixedGc;
    params.rho_q = rhoList(rIdx);
    params.calcGStat = false;
    params.edgeTrackPairs = singlePair;
    params.edgeTrackStep = recordStep;
    R = FHNres(x0,params);
    singleEdgeTime{rIdx} = R.edgeTrackTime;
    singleEdgeG{rIdx} = R.edgeTrackG;
end
params = params0;
params.gc = fixedGc;
params.rho_q = targetRho;
params.calcGStat = false;
params.edgeTrackPairs = trackPairs;
params.edgeTrackStep = recordStep;
R = FHNres(x0,params);
Fig4Data = struct();
Fig4Data.Gc = Gc;
Fig4Data.rhoList = rhoList;
Fig4Data.fixedGc = fixedGc;
Fig4Data.targetRho = targetRho;
Fig4Data.gEdgeCV = gEdgeCV;
Fig4Data.singleEdgePair = singlePair;
Fig4Data.singleEdgeTime = singleEdgeTime;
Fig4Data.singleEdgeG = singleEdgeG;
Fig4Data.trackEdgePairs = trackPairs;
Fig4Data.multiEdgeTime = R.edgeTrackTime;
Fig4Data.multiEdgeG = R.edgeTrackG;
save('Fig4.mat','Fig4Data','-v7.3')
toc

end
