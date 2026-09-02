function Fig3
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
numRepeat = 10;
Gc = (0.01:0.01:0.50).';
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
params0.couplingType = 'constant';
params0.rho_q = 0;
params0.epsE = 1e-5;
params0.q0 = ones(net.M_edge,1);
params0.qMin = 1e-12;
params0.calcCSF = false;
params0.calcHSF = true;
params0.calcEn = true;
params0.calcECF = true;
params0.calcHJ = true;
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
ResultsFlat = NaN(numG*numRepeat,5);
parfor job = 1:numG*numRepeat
    [gIdx,rIdx] = ind2sub([numG,numRepeat],job);
    params = params0;
    params.gc = Gc(gIdx);
    R = FHNres(x0Cell{rIdx},params);
    ResultsFlat(job,:) = [R.hSF,R.ecf,mean(R.En,'omitnan'),NodeEnergyCV(R.En),R.HJ];
end
Results = reshape(ResultsFlat,numG,numRepeat,5);
EcfData = Results(:,:,2);
EcfPositive = EcfData(isfinite(EcfData) & EcfData > 0);
if ~isempty(EcfPositive)
    Results(:,:,2) = EcfData/max(EcfPositive);
end
Mean = squeeze(mean(Results,2,'omitnan'));
Std = squeeze(std(Results,0,2,'omitnan'));
Fig3Data = struct();
Fig3Data.Gc = Gc;
Fig3Data.metricNames = {'HSF','ECF','MeanEnergy','EnergyCV','HJ'};
Fig3Data.numRepeat = numRepeat;
Fig3Data.N = N;
Fig3Data.K = K;
Fig3Data.p = p;
Fig3Data.all = Results;
Fig3Data.mean = Mean;
Fig3Data.std = Std;
Fig3Data.HSF.mean = Mean(:,1);
Fig3Data.HSF.std = Std(:,1);
Fig3Data.ECF.mean = Mean(:,2);
Fig3Data.ECF.std = Std(:,2);
Fig3Data.MeanEnergy.mean = Mean(:,3);
Fig3Data.MeanEnergy.std = Std(:,3);
Fig3Data.EnergyCV.mean = Mean(:,4);
Fig3Data.EnergyCV.std = Std(:,4);
Fig3Data.HJ.mean = Mean(:,5);
Fig3Data.HJ.std = Std(:,5);
save('Fig3.mat','Fig3Data','-v7.3')
toc

end
