function Fig5
clc
close all
tic
addpath('NetworkAnalyse')

N = 200;
K = 10;
p = 0.05;
networkSeed = 1;
numRepeat = 10;

GcBase = (0.01:0.01:0.40).';
GcControlExtra = (0.41:0.01:0.50).';
GcEnergyExtra = (0.41:0.01:0.45).';

rhoList = [0.05 0.20 0.35];
groupTypes = {'constant','energyExp','energyExp','energyExp'};
groupRho = [0 rhoList];

Net = SWStructureN(N,K,p,networkSeed);
net = BuildNetwork(Net.Adj);

params0 = struct( ...
    'h',0.01, ...
    'T',200001, ...
    'a',0.3, ...
    'b',1.0, ...
    'c',0.1, ...
    'apha',0.1, ...
    'beta',0.77, ...
    'k1',1.0, ...
    'k2',1.8);

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
params0.rescaleQMean = false;
params0.epsE = 1e-5;

params0.calcCSF = false;
params0.calcEn = false;
params0.calcHJ = false;
params0.calcGStat = false;
params0.calcFig7Map = false;

params0.calcEnergySyncTime = true;
params0.energySyncThreshold = 1e-5;
params0.energySyncHoldWindow = 1000;
params0.energySyncStartStep = params0.transient;

rng(1,'twister')
x0Cell = cell(numRepeat,1);
for r = 1:numRepeat
    x0Cell{r} = 2*rand(N,3);
end

if isempty(gcp('nocreate'))
    parpool;
end

numG = numel(GcBase);
numGroup = numel(groupTypes);

paramsBase = params0;
paramsBase.calcHSF = true;
paramsBase.calcECF = true;

ResultsBase = NaN(numG*numRepeat*numGroup,4);

parfor job = 1:numG*numRepeat*numGroup
    [gIdx,rIdx,mIdx] = ind2sub([numG,numRepeat,numGroup],job);

    params = paramsBase;
    params.gc = GcBase(gIdx);
    params.couplingType = groupTypes{mIdx};
    params.rho_q = groupRho(mIdx);

    R = FHNres(x0Cell{rIdx},params);

    ResultsBase(job,:) = [ ...
        R.hSF, ...
        R.ecf, ...
        R.energySyncTime, ...
        double(R.energySyncReached)];
end

AllBase = reshape(ResultsBase,numG,numRepeat,numGroup,4);

for mIdx = 1:numGroup
    E = AllBase(:,:,mIdx,2);
    v = E(isfinite(E) & E>0);
    if ~isempty(v)
        AllBase(:,:,mIdx,2) = E/max(v);
    end
end

paramsSync = params0;
paramsSync.calcHSF = false;
paramsSync.calcECF = false;

numGC = numel(GcControlExtra);
ResultsControl = NaN(numGC*numRepeat,2);

parfor job = 1:numGC*numRepeat
    [gIdx,rIdx] = ind2sub([numGC,numRepeat],job);

    params = paramsSync;
    params.gc = GcControlExtra(gIdx);
    params.couplingType = 'constant';
    params.rho_q = 0;

    R = FHNres(x0Cell{rIdx},params);

    ResultsControl(job,:) = [ ...
        R.energySyncTime, ...
        double(R.energySyncReached)];
end

AllControl = reshape(ResultsControl,numGC,numRepeat,2);

numGE = numel(GcEnergyExtra);
numEnergy = numel(rhoList);
ResultsEnergy = NaN(numGE*numRepeat*numEnergy,2);

parfor job = 1:numGE*numRepeat*numEnergy
    [gIdx,rIdx,mIdx] = ind2sub([numGE,numRepeat,numEnergy],job);

    params = paramsSync;
    params.gc = GcEnergyExtra(gIdx);
    params.couplingType = 'energyExp';
    params.rho_q = rhoList(mIdx);

    R = FHNres(x0Cell{rIdx},params);

    ResultsEnergy(job,:) = [ ...
        R.energySyncTime, ...
        double(R.energySyncReached)];
end

AllEnergy = reshape(ResultsEnergy,numGE,numRepeat,numEnergy,2);

Fig5Data = assembleFig5( ...
    GcBase, ...
    GcControlExtra, ...
    GcEnergyExtra, ...
    rhoList, ...
    groupTypes, ...
    groupRho, ...
    numRepeat, ...
    AllBase, ...
    AllControl, ...
    AllEnergy);

Fig5Data.N = N;
Fig5Data.K = K;
Fig5Data.p = p;
Fig5Data.networkSeed = networkSeed;
Fig5Data.h = params0.h;
Fig5Data.T = params0.T;
Fig5Data.transient = params0.transient;
Fig5Data.energySyncThreshold = params0.energySyncThreshold;
Fig5Data.energySyncHoldWindow = params0.energySyncHoldWindow;
Fig5Data.rescaleQMean = false;

save('Fig5.mat','Fig5Data','-v7.3')

fprintf('\nFig5 completed.\n')
fprintf('Final data saved to Fig5.mat\n')
fprintf('Critical gc: Control %.2f, sigma=0.05 %.2f, sigma=0.20 %.2f, sigma=0.35 %.2f\n', ...
    Fig5Data.criticalGc(1),Fig5Data.criticalGc(2), ...
    Fig5Data.criticalGc(3),Fig5Data.criticalGc(4))

toc
end


function Fig5Data = assembleFig5( ...
    GcBase,GcControlExtra,GcEnergyExtra,rhoList, ...
    groupTypes,groupRho,numRepeat,AllBase,AllControl,AllEnergy)

numGroup = numel(groupTypes);

HSFmean = squeeze(mean(AllBase(:,:,:,1),2,'omitnan'));
HSFstd = squeeze(std(AllBase(:,:,:,1),0,2,'omitnan'));

ECFmean = squeeze(mean(AllBase(:,:,:,2),2,'omitnan'));
ECFstd = squeeze(std(AllBase(:,:,:,2),0,2,'omitnan'));

GcSync = unique([GcBase(:);GcControlExtra(:);GcEnergyExtra(:)]);
GcSync = sort(GcSync);

SyncAll = NaN(numel(GcSync),numRepeat,numGroup,2);

for i = 1:numel(GcBase)
    idx = find(abs(GcSync-GcBase(i))<1e-12,1);
    SyncAll(idx,:,:,1) = AllBase(i,:,:,3);
    SyncAll(idx,:,:,2) = AllBase(i,:,:,4);
end

for i = 1:numel(GcControlExtra)
    idx = find(abs(GcSync-GcControlExtra(i))<1e-12,1);
    SyncAll(idx,:,1,1) = AllControl(i,:,1);
    SyncAll(idx,:,1,2) = AllControl(i,:,2);
end

for i = 1:numel(GcEnergyExtra)
    idx = find(abs(GcSync-GcEnergyExtra(i))<1e-12,1);
    for mIdx = 1:3
        SyncAll(idx,:,mIdx+1,1) = AllEnergy(i,:,mIdx,1);
        SyncAll(idx,:,mIdx+1,2) = AllEnergy(i,:,mIdx,2);
    end
end

Tmean = squeeze(mean(SyncAll(:,:,:,1),2,'omitnan'));
Tstd = squeeze(std(SyncAll(:,:,:,1),0,2,'omitnan'));
Psync = squeeze(mean(SyncAll(:,:,:,2),2,'omitnan'));
nSync = squeeze(sum(SyncAll(:,:,:,2)>0.5,2));

criticalGc = NaN(numGroup,1);
for mIdx = 1:numGroup
    idx = find(Psync(:,mIdx)>=1,1,'first');
    if ~isempty(idx)
        criticalGc(mIdx) = GcSync(idx);
    end
end

Fig5Data = struct();

Fig5Data.Gc = GcBase(:);
Fig5Data.GcSync = GcSync(:);
Fig5Data.GcSyncDisplayMax = 0.45;

Fig5Data.rhoList = rhoList;
Fig5Data.groupTypes = groupTypes;
Fig5Data.groupRho = groupRho;
Fig5Data.numRepeat = numRepeat;

Fig5Data.HSF.mean = HSFmean;
Fig5Data.HSF.std = HSFstd;

Fig5Data.ECF.mean = ECFmean;
Fig5Data.ECF.std = ECFstd;

Fig5Data.Tsyn.mean = Tmean;
Fig5Data.Tsyn.std = Tstd;

Fig5Data.Psync = Psync;
Fig5Data.nSync = nSync;
Fig5Data.criticalGc = criticalGc;

Fig5Data.all = AllBase;
Fig5Data.syncAll = SyncAll;

Fig5Data.range.base = [min(GcBase) max(GcBase)];
Fig5Data.range.controlExtra = [min(GcControlExtra) max(GcControlExtra)];
Fig5Data.range.energyExtra = [min(GcEnergyExtra) max(GcEnergyExtra)];
end
