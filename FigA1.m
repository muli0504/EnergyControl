function FigA1
clc
close all
tic
addpath('NetworkAnalyse')
N = 200;
K = 10;
p = 0.05;
numRepeat = 10;
Gc = (0.01:0.01:0.40).';
rho = 0.35;
gcDetail = 0.30;
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
params0.couplingType = 'energyExp';
params0.rho_q = rho;
params0.epsE = 1e-5;
params0.rescaleQMean = false;
params0.calcCSF = true;
params0.calcHSF = true;
params0.calcEn = false;
params0.calcECF = false;
params0.calcHJ = false;
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
Results = NaN(numG*numRepeat,2);
parfor job = 1:numG*numRepeat
    [gIdx,rIdx] = ind2sub([numG,numRepeat],job);
    params = params0;
    params.gc = Gc(gIdx);
    R = FHNres(x0Cell{rIdx},params);
    Results(job,:) = [R.hSF,R.CSF];
end
All = reshape(Results,numG,numRepeat,2);
[~,idxDetail] = min(abs(Gc-gcDetail));
paramsDetail = params0;
paramsDetail.gc = Gc(idxDetail);
paramsDetail.detailRecordStep = 10;
Detail = FHNdetailErrors_RMS(x0Cell{1},paramsDetail);
FigA1Data = struct();
FigA1Data.Gc = Gc;
FigA1Data.HSF.mean = mean(All(:,:,1),2,'omitnan');
FigA1Data.HSF.std = std(All(:,:,1),0,2,'omitnan');
FigA1Data.CSF.mean = mean(All(:,:,2),2,'omitnan');
FigA1Data.CSF.std = std(All(:,:,2),0,2,'omitnan');
FigA1Data.gcDetail = Gc(idxDetail);
FigA1Data.Detail = Detail;
FigA1Data.all = All;
save('FigA1.mat','FigA1Data','-v7.3')
toc

end
