function Fig9
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
numRepeat = 10;
Gc = (0.001:0.001:0.10).';
rhoList = [0.05 0.2 0.35];
groupTypes = {'constant','energyExp','energyExp','energyExp'};
groupRho = [0 rhoList];
params0 = struct('h',0.01,'T',200001,'a',0.3,'b',1.0,'c',0.1,'apha',0.1,'beta',0.77,'k1',1.0,'k2',1.8);
params0.Amp = 0.828*ones(N,1);
params0.freq = 0.19*ones(N,1);
params0.transient = round(0.3*params0.T);
params0.qMin = 1e-12;
params0.epsE = 1e-5;
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
if isempty(gcp('nocreate'))
    parpool;
end
SF = runTopology(SFStructureN(N),params0,Gc,groupTypes,groupRho,numRepeat);
ER = runTopology(ERStructureN(N),params0,Gc,groupTypes,groupRho,numRepeat);
Fig9Data = struct();
Fig9Data.Gc = Gc;
Fig9Data.rhoList = rhoList;
Fig9Data.groupTypes = groupTypes;
Fig9Data.groupRho = groupRho;
Fig9Data.SF = SF;
Fig9Data.ER = ER;
save('Fig9.mat','Fig9Data','-v7.3')
toc


end
function S = runTopology(Net,params0,Gc,groupTypes,groupRho,numRepeat)
net = BuildNetwork(Net.Adj);
N = net.N;
params0.N = N;
params0.Adj = net.Adj;
params0.deg = net.deg;
params0.edge_i = net.edge_i;
params0.edge_j = net.edge_j;
params0.M_edge = net.M_edge;
params0.B = net.B;
params0.q0 = ones(net.M_edge,1);
rng(1,'twister')
x0Cell = cell(numRepeat,1);
for r = 1:numRepeat
    x0Cell{r} = 2*rand(N,3);
end
numG = numel(Gc);
numGroup = numel(groupTypes);
Results = NaN(numG*numRepeat*numGroup,4);
parfor job = 1:numG*numRepeat*numGroup
    [gIdx,rIdx,mIdx] = ind2sub([numG,numRepeat,numGroup],job);
    params = params0;
    params.gc = Gc(gIdx);
    params.couplingType = groupTypes{mIdx};
    params.rho_q = groupRho(mIdx);
    R = FHNres(x0Cell{rIdx},params);
    Results(job,:) = [R.hSF,R.energySyncTime,R.HJ,double(R.energySyncReached)];
end
All = reshape(Results,numG,numRepeat,numGroup,4);
S.HSF.mean = squeeze(mean(All(:,:,:,1),2,'omitnan'));
S.HSF.std = squeeze(std(All(:,:,:,1),0,2,'omitnan'));
S.Tsyn.mean = squeeze(mean(All(:,:,:,2),2,'omitnan'));
S.Tsyn.std = squeeze(std(All(:,:,:,2),0,2,'omitnan'));
S.HJ.mean = squeeze(mean(All(:,:,:,3),2,'omitnan'));
S.HJ.std = squeeze(std(All(:,:,:,3),0,2,'omitnan'));
S.Psync = squeeze(mean(All(:,:,:,4),2,'omitnan'));
S.criticalGc = NaN(numGroup,1);
for mIdx = 1:numGroup
    idx = find(S.Psync(:,mIdx)>=1,1,'first');
    if ~isempty(idx)
        S.criticalGc(mIdx) = Gc(idx);
    end
end
S.all = All;
end
