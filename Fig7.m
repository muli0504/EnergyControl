function Fig7
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
fixedGc = 0.22;
rhoList = [0.05 0.2 0.35];
mapStep = 10;
Net = SWStructureN(N,K,p,1);
net = BuildNetwork(Net.Adj);
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
params0.gc = fixedGc;
params0.couplingType = 'energyExp';
params0.epsE = 1e-5;
params0.calcCSF = false;
params0.calcHSF = false;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = false;
params0.calcGStat = false;
params0.calcEnergySyncTime = false;
params0.calcFig7Map = true;
params0.fig7MapStep = mapStep;
rng(1,'twister')
x0 = 2*rand(N,3);
Fig7Data = struct();
Fig7Data.gc = fixedGc;
Fig7Data.rhoList = rhoList;
Fig7Data.edgePairs = [net.edge_i,net.edge_j];
Fig7Data.time = cell(numel(rhoList),1);
Fig7Data.VmMap = cell(numel(rhoList),1);
Fig7Data.EnergyMap = cell(numel(rhoList),1);
Fig7Data.edgeHJMap = cell(numel(rhoList),1);
Fig7Data.edgeGMap = cell(numel(rhoList),1);
for rIdx = 1:numel(rhoList)
    params = params0;
    params.rho_q = rhoList(rIdx);
    R = FHNres(x0,params);
    Fig7Data.time{rIdx} = R.mapTime;
    Fig7Data.VmMap{rIdx} = R.VmMap;
    Fig7Data.EnergyMap{rIdx} = R.EnergyMap;
    Fig7Data.edgeHJMap{rIdx} = R.edgeHJMap;
    Fig7Data.edgeGMap{rIdx} = R.edgeGMap;
end
save('Fig7.mat','Fig7Data','-v7.3')
toc

end
