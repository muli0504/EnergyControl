function FigA3
clc
close all
tic
addpath('NetworkAnalyse')
numRepeat = 10;
gc = 0.30;
rho = 0.35;
epsE = 1e-5;
NList = [200 500 1000 2000 3000 4000];
Kfixed = 10;
p = 0.05;
NResult = runN(NList,Kfixed,p,numRepeat,gc,rho,epsE);
KList = [4 6 8 10 12 14 16 18 20];
Nfixed = 200;
KResult = runK(KList,Nfixed,p,numRepeat,gc,rho,epsE);
FigA3Data = struct();
FigA3Data.N = NResult;
FigA3Data.K = KResult;
FigA3Data.gc = gc;
FigA3Data.rho = rho;
FigA3Data.epsE = epsE;
save('FigA3.mat','FigA3Data','-v7.3')
toc


end
function S = runN(NList,K,p,numRepeat,gc,rho,epsE)
numN = numel(NList);
HSF = NaN(numN,numRepeat,2);
HJ = NaN(numN,numRepeat,2);
edges = NaN(numN,1);
for nIdx = 1:numN
    N = NList(nIdx);
    Net = SWStructureN(N,K,p,1);
    net = BuildNetwork(Net.Adj);
    edges(nIdx) = net.M_edge;
    params0 = baseParams(N,net,gc,epsE);
    rng(1000+nIdx,'twister')
    x0Cell = cell(numRepeat,1);
    for r = 1:numRepeat
        x0Cell{r} = 2*rand(N,3);
    end
    Results = NaN(numRepeat*2,2);
    parfor job = 1:numRepeat*2
        [rIdx,mIdx] = ind2sub([numRepeat,2],job);
        params = params0;
        if mIdx == 1
            params.couplingType = 'constant';
            params.rho_q = 0;
        else
            params.couplingType = 'energyExp';
            params.rho_q = rho;
        end
        R = FHNres(x0Cell{rIdx},params);
        Results(job,:) = [R.hSF,R.HJ];
    end
    R3 = reshape(Results,numRepeat,2,2);
    HSF(nIdx,:,:) = permute(R3(:,:,1),[3 1 2]);
    HJ(nIdx,:,:) = permute(R3(:,:,2),[3 1 2]);
end
S.values = NList(:);
S.edges = edges;
S.HSF.mean = squeeze(mean(HSF,2,'omitnan'));
S.HSF.std = squeeze(std(HSF,0,2,'omitnan'));
S.HJ.mean = squeeze(mean(HJ,2,'omitnan'));
S.HJ.std = squeeze(std(HJ,0,2,'omitnan'));
S.HSF.all = HSF;
S.HJ.all = HJ;
end

function S = runK(KList,N,p,numRepeat,gc,rho,epsE)
numK = numel(KList);
HSF = NaN(numK,numRepeat,2);
HJ = NaN(numK,numRepeat,2);
edges = NaN(numK,1);
for kIdx = 1:numK
    K = KList(kIdx);
    Net = SWStructureN(N,K,p,1);
    net = BuildNetwork(Net.Adj);
    edges(kIdx) = net.M_edge;
    params0 = baseParams(N,net,gc,epsE);
    rng(1000+kIdx,'twister')
    x0Cell = cell(numRepeat,1);
    for r = 1:numRepeat
        x0Cell{r} = 2*rand(N,3);
    end
    Results = NaN(numRepeat*2,2);
    parfor job = 1:numRepeat*2
        [rIdx,mIdx] = ind2sub([numRepeat,2],job);
        params = params0;
        if mIdx == 1
            params.couplingType = 'constant';
            params.rho_q = 0;
        else
            params.couplingType = 'energyExp';
            params.rho_q = rho;
        end
        R = FHNres(x0Cell{rIdx},params);
        Results(job,:) = [R.hSF,R.HJ];
    end
    R3 = reshape(Results,numRepeat,2,2);
    HSF(kIdx,:,:) = permute(R3(:,:,1),[3 1 2]);
    HJ(kIdx,:,:) = permute(R3(:,:,2),[3 1 2]);
end
S.values = KList(:);
S.edges = edges;
S.HSF.mean = squeeze(mean(HSF,2,'omitnan'));
S.HSF.std = squeeze(std(HSF,0,2,'omitnan'));
S.HJ.mean = squeeze(mean(HJ,2,'omitnan'));
S.HJ.std = squeeze(std(HJ,0,2,'omitnan'));
S.HSF.all = HSF;
S.HJ.all = HJ;
end

function params = baseParams(N,net,gc,epsE)
params = struct('h',0.01,'T',200001,'a',0.3,'b',1.0,'c',0.1,'apha',0.1,'beta',0.77,'k1',1.0,'k2',1.8);
params.N = N;
params.Adj = net.Adj;
params.deg = net.deg;
params.edge_i = net.edge_i;
params.edge_j = net.edge_j;
params.M_edge = net.M_edge;
params.B = net.B;
params.Amp = 0.828*ones(N,1);
params.freq = 0.19*ones(N,1);
params.transient = round(0.3*params.T);
params.q0 = ones(net.M_edge,1);
params.qMin = 1e-12;
params.gc = gc;
params.epsE = epsE;
params.calcCSF = false;
params.calcHSF = true;
params.calcEn = false;
params.calcECF = false;
params.calcHJ = true;
params.calcGStat = false;
params.calcEnergySyncTime = false;
params.calcFig7Map = false;
end
