function Fig11
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
numRepeat = 10;
Gc = (0.01:0.01:0.45).';
rho = 0.35;
groupTypes = {'energyExp','hebbianPRX','pddp'};
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
params0.qBase = 1;
params0.pHebb = 2.5;
params0.kHebb = rho;
params0.rho_pddp = rho;
params0.pddpAmp = 0.90;
[params0.phaseCenterX,params0.phaseCenterY] = PhaseCenter(params0,0.828,0.19,200000,0.3);
params0.calcCSF = false;
params0.calcHSF = true;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = true;
params0.calcGStat = false;
params0.calcEnergySyncTime = true;
params0.energySyncThreshold = 1e-5;
params0.energySyncHoldWindow = 1000;
params0.energySyncStartStep = params0.transient;
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
Results = NaN(numG*numRepeat*numGroup,4);
parfor job = 1:numG*numRepeat*numGroup
    [gIdx,rIdx,mIdx] = ind2sub([numG,numRepeat,numGroup],job);
    params = params0;
    params.gc = Gc(gIdx);
    params.couplingType = groupTypes{mIdx};
    params.rho_q = rho;
    R = FHNres(x0Cell{rIdx},params);
    Results(job,:) = [R.hSF,R.HJ,R.energySyncTime,double(R.energySyncReached)];
end
All = reshape(Results,numG,numRepeat,numGroup,4);
HSFmean = squeeze(mean(All(:,:,:,1),2,'omitnan'));
HSFstd = squeeze(std(All(:,:,:,1),0,2,'omitnan'));
HJmean = squeeze(mean(All(:,:,:,2),2,'omitnan'));
HJstd = squeeze(std(All(:,:,:,2),0,2,'omitnan'));
Tmean = squeeze(mean(All(:,:,:,3),2,'omitnan'));
Tstd = squeeze(std(All(:,:,:,3),0,2,'omitnan'));
Psync = squeeze(mean(All(:,:,:,4),2,'omitnan'));
criticalGc = NaN(numGroup,1);
for mIdx = 1:numGroup
    idx = find(Psync(:,mIdx)>=1,1,'first');
    if ~isempty(idx)
        criticalGc(mIdx) = Gc(idx);
    end
end
Fig11Data = struct();
Fig11Data.Gc = Gc;
Fig11Data.groupTypes = groupTypes;
Fig11Data.rho = rho;
Fig11Data.HSF.mean = HSFmean;
Fig11Data.HSF.std = HSFstd;
Fig11Data.HJ.mean = HJmean;
Fig11Data.HJ.std = HJstd;
Fig11Data.Tsyn.mean = Tmean;
Fig11Data.Tsyn.std = Tstd;
Fig11Data.Psync = Psync;
Fig11Data.criticalGc = criticalGc;
Fig11Data.phaseCenter = [params0.phaseCenterX params0.phaseCenterY];
Fig11Data.all = All;
save('Fig11.mat','Fig11Data','-v7.3')
toc

end
