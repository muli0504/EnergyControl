function Fig8
clc
close all
tic
addpath('NetworkAnalyse')

N = 200;
K = 10;
p = 0.05;
networkSeed = 1;
numRepeat = 10;

Gc = (0.01:0.01:0.40).';
gcCompare = 0.30;
[~,gcCompareIdx] = min(abs(Gc-gcCompare));

rho = 0.35;
epsE = 1e-5;

betaList = [0.73 0.74 0.75 0.77 0.85 0.87 0.89 0.92].';
fList = [0.179 0.181 0.190 0.192 0.193 0.194 0.195 0.196].';

Net = SWStructureN(N,K,p,networkSeed);
net = BuildNetwork(Net.Adj);

params0 = struct();
params0.h = 0.01;
params0.T = 200001;
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
params0.transient = round(0.3*params0.T);
params0.q0 = ones(net.M_edge,1);
params0.qMin = 1e-12;
params0.epsE = epsE;
params0.calcCSF = false;
params0.calcHSF = true;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = true;
params0.calcGStat = false;
params0.calcEnergySyncTime = true;
params0.energySyncThreshold = epsE;
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

Beta = runSensitivity(betaList,'beta',rho,0.19,Gc,gcCompareIdx,x0Cell,params0);
Freq = runSensitivity(fList,'freq',rho,0.77,Gc,gcCompareIdx,x0Cell,params0);

Fig8Data = struct();
Fig8Data.Gc = Gc;
Fig8Data.gcCompare = Gc(gcCompareIdx);
Fig8Data.rho = rho;
Fig8Data.epsE = epsE;
Fig8Data.numRepeat = numRepeat;
Fig8Data.beta = Beta;
Fig8Data.frequency = Freq;

save('Fig8.mat','Fig8Data','-v7.3')
toc


end
function S = runSensitivity(values,kind,rho,fixedValue,Gc,gcCompareIdx,x0Cell,params0)

numV = numel(values);
numG = numel(Gc);
numRepeat = numel(x0Cell);
groupTypes = {'constant','energyExp'};
groupRho = [0 rho];
numGroup = 2;
All = NaN(numV,numG,numRepeat,numGroup,4);

for vIdx = 1:numV
    Results = NaN(numG*numRepeat*numGroup,4);

    parfor job = 1:numG*numRepeat*numGroup
        [gIdx,rIdx,mIdx] = ind2sub([numG,numRepeat,numGroup],job);

        params = params0;
        params.gc = Gc(gIdx);
        params.couplingType = groupTypes{mIdx};
        params.rho_q = groupRho(mIdx);

        if strcmp(kind,'beta')
            params.beta = values(vIdx);
            params.freq = fixedValue*ones(params.N,1);
        else
            params.beta = fixedValue;
            params.freq = values(vIdx)*ones(params.N,1);
        end

        R = FHNres(x0Cell{rIdx},params);

        Results(job,:) = [ ...
            R.hSF, ...
            R.energySyncTime, ...
            R.HJ, ...
            double(R.energySyncReached)];
    end

    All(vIdx,:,:,:,:) = reshape(Results,1,numG,numRepeat,numGroup,4);
end

HSFmean = squeeze(mean(All(:,:,:,:,1),3,'omitnan'));
HSFstd = squeeze(std(All(:,:,:,:,1),0,3,'omitnan'));

TsynMean = squeeze(mean(All(:,:,:,:,2),3,'omitnan'));
TsynStd = squeeze(std(All(:,:,:,:,2),0,3,'omitnan'));

HJmean = squeeze(mean(All(:,:,:,:,3),3,'omitnan'));
HJstd = squeeze(std(All(:,:,:,:,3),0,3,'omitnan'));

Psync = squeeze(mean(All(:,:,:,:,4),3,'omitnan'));

criticalGc = NaN(numV,numGroup);
for vIdx = 1:numV
    for mIdx = 1:numGroup
        idx = find(Psync(vIdx,:,mIdx)>=1,1,'first');
        if ~isempty(idx)
            criticalGc(vIdx,mIdx) = Gc(idx);
        end
    end
end

Tmax = params0.T*params0.h;
TsynPlot = TsynMean;
TsynPlot(~isfinite(TsynPlot)) = Tmax;

S.values = values;
S.rho = rho;
S.groupTypes = groupTypes;
S.criticalGc = criticalGc;
S.HSF_gc030 = squeeze(HSFmean(:,gcCompareIdx,:));
S.HSF_std_gc030 = squeeze(HSFstd(:,gcCompareIdx,:));
S.Tsyn_gc030 = squeeze(TsynMean(:,gcCompareIdx,:));
S.Tsyn_std_gc030 = squeeze(TsynStd(:,gcCompareIdx,:));
S.TsynPlot_gc030 = squeeze(TsynPlot(:,gcCompareIdx,:));
S.HJ_gc030 = squeeze(HJmean(:,gcCompareIdx,:));
S.HJ_std_gc030 = squeeze(HJstd(:,gcCompareIdx,:));
S.Psync = Psync;
S.all = All;

end
