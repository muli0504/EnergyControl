function TableRobustness
clc
close all
tic

addpath('NetworkAnalyse')

N = 200;
K = 10;
p = 0.05;
networkSeed = 1;
numRepeat = 10;

Gc = (0.01:0.01:0.45).';
rho = 0.35;
gcCompare = 0.30;

Net = SWStructureN(N,K,p,networkSeed);
net = BuildNetwork(Net.Adj);

params0 = struct();
params0.a = 0.3;
params0.b = 1.0;
params0.c = 0.1;
params0.apha = 0.1;
params0.beta = 0.77;
params0.k1 = 1.0;
params0.k2 = 1.8;

params0.N = N;
params0.Adj = net.Adj;
params0.deg = net.deg;
params0.edge_i = net.edge_i;
params0.edge_j = net.edge_j;
params0.M_edge = net.M_edge;
params0.B = net.B;

params0.Amp = 0.828*ones(N,1);
params0.freq = 0.19*ones(N,1);

params0.q0 = ones(net.M_edge,1);
params0.qMin = 1e-12;

params0.calcCSF = false;
params0.calcHSF = false;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = true;
params0.calcGStat = false;
params0.calcEnergySyncTime = true;
params0.calcFig7Map = false;

rng(1,'twister')

x0Cell = cell(numRepeat,1);
for r = 1:numRepeat
    x0Cell{r} = 2*rand(N,3);
end

if isempty(gcp('nocreate'))
    parpool;
end

H = runH(params0,x0Cell,Gc,rho,gcCompare);
EpsH = runEpsH(params0,x0Cell,Gc,rho,gcCompare);
Tsim = runTsim(params0,x0Cell,Gc,rho,gcCompare);
Q0 = runQ0(params0,x0Cell,Gc,rho,gcCompare);

TableRobustnessData = struct();
TableRobustnessData.N = N;
TableRobustnessData.K = K;
TableRobustnessData.p = p;
TableRobustnessData.networkSeed = networkSeed;
TableRobustnessData.numRepeat = numRepeat;
TableRobustnessData.Gc = Gc;
TableRobustnessData.rho = rho;
TableRobustnessData.gcCompare = gcCompare;
TableRobustnessData.H = H;
TableRobustnessData.EpsH = EpsH;
TableRobustnessData.Tsim = Tsim;
TableRobustnessData.Q0 = Q0;

save('TableRobustness.mat','TableRobustnessData','-v7.3')

toc
end

function S = runH(params0,x0Cell,Gc,rho,gcCompare)

values = [0.02 0.01 0.005];
numC = numel(values);
numG = numel(Gc);
numR = numel(x0Cell);
numM = 2;

Rall = NaN(numG,numR,numM,numC,3);

for cIdx = 1:numC

    h = values(cIdx);
    Results = NaN(numG*numR*numM,3);

    parfor job = 1:numG*numR*numM

        [gIdx,rIdx,mIdx] = ind2sub([numG,numR,numM],job);

        params = params0;
        params.h = h;
        params.T = round(2000/h)+1;
        params.transient = round(600/h);
        params.energySyncStartStep = params.transient;
        params.energySyncThreshold = 1e-5;
        params.energySyncHoldWindow = round(10/h);
        params.epsE = 1e-5;
        params.gc = Gc(gIdx);

        if mIdx == 1
            params.couplingType = 'constant';
            params.rho_q = 0;
        else
            params.couplingType = 'energyExp';
            params.rho_q = rho;
        end

        R = FHNres(x0Cell{rIdx},params);

        Results(job,:) = [ ...
            R.energySyncTime, ...
            double(R.energySyncReached), ...
            R.HJ];
    end

    Rall(:,:,:,cIdx,:) = reshape(Results,numG,numR,numM,1,3);
end

S = reduceGrid(Rall,Gc,gcCompare,1.0);
S.values = values;
S.Tsim = 2000;
S.transientTime = 600;
S.holdTime = 10;
S.epsH = 1e-5;
S.epsSyn = 1e-5;
S.all = Rall;

end

function S = runEpsH(params0,x0Cell,Gc,rho,gcCompare)

values = [1e-6 5e-6 1e-5];
numC = numel(values);
numG = numel(Gc);
numR = numel(x0Cell);
numM = 2;

Rall = NaN(numG,numR,numM,numC,3);

for cIdx = 1:numC

    Results = NaN(numG*numR*numM,3);

    parfor job = 1:numG*numR*numM

        [gIdx,rIdx,mIdx] = ind2sub([numG,numR,numM],job);

        params = params0;
        params.h = 0.01;
        params.T = 200001;
        params.transient = round(0.3*params.T);
        params.energySyncStartStep = params.transient;
        params.energySyncThreshold = 1e-5;
        params.energySyncHoldWindow = 1000;
        params.epsE = values(cIdx);
        params.gc = Gc(gIdx);

        if mIdx == 1
            params.couplingType = 'constant';
            params.rho_q = 0;
        else
            params.couplingType = 'energyExp';
            params.rho_q = rho;
        end

        R = FHNres(x0Cell{rIdx},params);

        Results(job,:) = [ ...
            R.energySyncTime, ...
            double(R.energySyncReached), ...
            R.HJ];
    end

    Rall(:,:,:,cIdx,:) = reshape(Results,numG,numR,numM,1,3);
end

S = reduceGrid(Rall,Gc,gcCompare,1.0);
S.values = values;
S.epsSyn = 1e-5;
S.h = 0.01;
S.T = 200001;
S.all = Rall;

end

function S = runTsim(params0,x0Cell,Gc,rho,gcCompare)

values = [1500 2000 2300];
numC = numel(values);
numG = numel(Gc);
numR = numel(x0Cell);
numM = 2;

Rall = NaN(numG,numR,numM,numC,3);

for cIdx = 1:numC

    Results = NaN(numG*numR*numM,3);

    parfor job = 1:numG*numR*numM

        [gIdx,rIdx,mIdx] = ind2sub([numG,numR,numM],job);

        params = params0;
        params.h = 0.01;
        params.T = round(values(cIdx)/params.h)+1;
        params.transient = round(600/params.h);
        params.energySyncStartStep = params.transient;
        params.energySyncThreshold = 1e-5;
        params.energySyncHoldWindow = round(10/params.h);
        params.epsE = 1e-5;
        params.rescaleQMean = true;
        params.gc = Gc(gIdx);

        if mIdx == 1
            params.couplingType = 'constant';
            params.rho_q = 0;
        else
            params.couplingType = 'energyExp';
            params.rho_q = rho;
        end

        R = FHNres(x0Cell{rIdx},params);

        Results(job,:) = [ ...
            R.energySyncTime, ...
            double(R.energySyncReached), ...
            R.HJ];
    end

    Rall(:,:,:,cIdx,:) = reshape(Results,numG,numR,numM,1,3);
end

S = reduceGrid(Rall,Gc,gcCompare,1.0);
S.values = values;
S.h = 0.01;
S.transientTime = 600;
S.holdTime = 10;
S.epsH = 1e-5;
S.epsSyn = 1e-5;
S.rescaleQMean = true;
S.all = Rall;

end

function S = runQ0(params0,x0Cell,Gc,rho,gcCompare)

q0List = [0.5 1.0 2.0];
numQ0 = numel(q0List);
numRepeat = numel(x0Cell);
nGc = numel(Gc);
numJobs = numQ0*numRepeat*nGc;

params0.h = 0.01;
Tsim = 2000;
params0.T = round(Tsim/params0.h);
params0.transient = round(600/params0.h);
params0.energySyncStartStep = params0.transient;
params0.energySyncThreshold = 1e-5;
params0.energySyncHoldWindow = round(10/params0.h);
params0.epsE = 1e-5;
params0.rescaleQMean = true;
params0.couplingType = 'energyExp';
params0.rho_q = rho;

Results = NaN(numJobs,3);

parfor job = 1:numJobs

    [qIdx,repIdx,gIdx] = ind2sub([numQ0,numRepeat,nGc],job);

    params = params0;
    params.gc = Gc(gIdx);
    params.q0 = q0List(qIdx)*ones(params.M_edge,1);

    R = FHNres(x0Cell{repIdx},params);

    Results(job,:) = [ ...
        R.energySyncTime, ...
        double(R.energySyncReached), ...
        R.HJ];
end

AllRaw = reshape(Results,numQ0,numRepeat,nGc,3);
All = permute(AllRaw,[3 2 1 4]);

TsyncMean = squeeze(mean(All(:,:,:,1),2,'omitnan'));
TsyncStd = squeeze(std(All(:,:,:,1),0,2,'omitnan'));

SyncRatio = squeeze(mean(All(:,:,:,2),2,'omitnan'));
nSync = squeeze(sum(All(:,:,:,2)>0.5,2));

HJmean = squeeze(mean(All(:,:,:,3),2,'omitnan'));
HJstd = squeeze(std(All(:,:,:,3),0,2,'omitnan'));

P0 = 0.90;
criticalGc = NaN(numQ0,1);

for qIdx = 1:numQ0
    idx = find(SyncRatio(:,qIdx)>=P0,1,'first');

    if ~isempty(idx)
        criticalGc(qIdx) = Gc(idx);
    end
end

[~,idxCompare] = min(abs(Gc-gcCompare));

S.values = q0List;
S.P0 = P0;
S.controlApplicable = false;
S.criticalGc = criticalGc;

S.Tsyn_gc030 = TsyncMean(idxCompare,:).';
S.Tsyn_std_gc030 = TsyncStd(idxCompare,:).';

S.HJ_gc030 = HJmean(idxCompare,:).';
S.HJ_std_gc030 = HJstd(idxCompare,:).';

S.SyncRatio_gc030 = SyncRatio(idxCompare,:).';
S.nSync_gc030 = nSync(idxCompare,:).';

S.SyncRatio = SyncRatio;
S.TsynMean = TsyncMean;
S.TsynStd = TsyncStd;
S.HJMean = HJmean;
S.HJStd = HJstd;

S.h = params0.h;
S.Tsim = Tsim;
S.transientTime = 600;
S.holdTime = 10;
S.epsH = 1e-5;
S.rescaleQMean = true;

S.ResultsRaw = Results;
S.All = All;

end

function S = reduceGrid(All,Gc,gcCompare,criticalRatio)

syncTimeMean = squeeze(mean(All(:,:,:,:,1),2,'omitnan'));
syncTimeStd = squeeze(std(All(:,:,:,:,1),0,2,'omitnan'));

syncRatio = squeeze(mean(All(:,:,:,:,2),2,'omitnan'));
nSync = squeeze(sum(All(:,:,:,:,2)>0.5,2));

HJMean = squeeze(mean(All(:,:,:,:,3),2,'omitnan'));
HJStd = squeeze(std(All(:,:,:,:,3),0,2,'omitnan'));

numM = size(syncRatio,2);
numC = size(syncRatio,3);

criticalGc = NaN(numM,numC);

for cIdx = 1:numC
    for mIdx = 1:numM

        idx = find(syncRatio(:,mIdx,cIdx)>=criticalRatio,1,'first');

        if ~isempty(idx)
            criticalGc(mIdx,cIdx) = Gc(idx);
        end
    end
end

[~,idxCompare] = min(abs(Gc-gcCompare));

S.criticalRatio = criticalRatio;
S.criticalGc = criticalGc;

S.Tsyn_gc030 = squeeze(syncTimeMean(idxCompare,:,:));
S.Tsyn_std_gc030 = squeeze(syncTimeStd(idxCompare,:,:));

S.HJ_gc030 = squeeze(HJMean(idxCompare,:,:));
S.HJ_std_gc030 = squeeze(HJStd(idxCompare,:,:));

S.SyncRatio_gc030 = squeeze(syncRatio(idxCompare,:,:));
S.nSync_gc030 = squeeze(nSync(idxCompare,:,:));

S.SyncRatio = syncRatio;
S.nSync = nSync;
S.TsynMean = syncTimeMean;
S.TsynStd = syncTimeStd;
S.HJMean = HJMean;
S.HJStd = HJStd;

end
