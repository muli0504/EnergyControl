function Fig6
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
numRepeat = 10;
Gc = (0.01:0.01:0.40).';
rhoList = [0.05 0.2 0.35];
groupTypes = {'constant','energyExp','energyExp','energyExp'};
groupRho = [0 rhoList];
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
params0.epsE = 1e-5;
params0.calcCSF = false;
params0.calcHSF = false;
params0.calcEn = true;
params0.calcECF = false;
params0.calcHJ = true;
params0.calcGStat = false;
params0.calcEnergySyncTime = false;
params0.calcFig7Map = false;
rng(1,'twister')
x0Cell = cell(numRepeat,1);
for r = 1:numRepeat
    x0Cell{r} = 2*rand(N,3);
end
if isempty(gcp('nocreate'))
    parpool;
end
numG = numel(Gc);
numGroup = numel(groupTypes);
ResultsFlat = NaN(numG*numRepeat*numGroup,3);
parfor job = 1:numG*numRepeat*numGroup
    [gIdx,rIdx,mIdx] = ind2sub([numG,numRepeat,numGroup],job);
    params = params0;
    params.gc = Gc(gIdx);
    params.couplingType = groupTypes{mIdx};
    params.rho_q = groupRho(mIdx);
    R = FHNres(x0Cell{rIdx},params);
    ResultsFlat(job,:) = [mean(R.En,'omitnan'),NodeEnergyCV(R.En),R.HJ];
end
All = reshape(ResultsFlat,numG,numRepeat,numGroup,3);
Mean = squeeze(mean(All,2,'omitnan'));
Std = squeeze(std(All,0,2,'omitnan'));
Fig6Data = struct();
Fig6Data.Gc = Gc;
Fig6Data.rhoList = rhoList;
Fig6Data.groupTypes = groupTypes;
Fig6Data.groupRho = groupRho;
Fig6Data.numRepeat = numRepeat;
Fig6Data.MeanEnergy.mean = Mean(:,:,1);
Fig6Data.MeanEnergy.std = Std(:,:,1);
Fig6Data.EnergyCV.mean = Mean(:,:,2);
Fig6Data.EnergyCV.std = Std(:,:,2);
Fig6Data.HJ.mean = Mean(:,:,3);
Fig6Data.HJ.std = Std(:,:,3);
Fig6Data.all = All;
save('Fig6.mat','Fig6Data','-v7.3')
toc

end
