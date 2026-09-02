function Fig10
clc
close all
tic
addpath('NetworkAnalyse')

N = 200;
K = 10;
p = 0.05;
networkSeed = 1;
numRepeat = 10;

GcPresent = (0.005:0.005:0.38).';
GcRef54 = (0.005:0.005:0.14).';

sigmaPresent = 0.001;
sigmaRef54 = 0.001;
epsH = 1e-5;

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
params0.epsH = epsH;
params0.transient = round(0.3*params0.T);
params0.qMin = 1e-12;
params0.kMin = 1e-12;
params0.maxAbsState = 1e8;
params0.maxAdaptive = 1e10;
params0.sigmaPresent = sigmaPresent;
params0.sigmaRef54 = sigmaRef54;

Net = SWStructureN(N,K,p,networkSeed);
net = BuildNetwork(Net.Adj);

params0.N = N;
params0.Adj = net.Adj;
params0.deg = net.deg;
params0.edge_i = net.edge_i;
params0.edge_j = net.edge_j;
params0.B = net.B;
params0.M_edge = net.M_edge;
params0.Amp = 0.828*ones(N,1);
params0.freq = 0.19*ones(N,1);

rng(1,'twister')
x0Cell = cell(numRepeat,1);
for r = 1:numRepeat
    x0Cell{r} = 2*rand(N,3);
end

if isempty(gcp('nocreate'))
    parpool;
end

Present = runMethod('present',GcPresent,x0Cell,params0);
Ref54 = runMethod('ref54',GcRef54,x0Cell,params0);

Fig10Data = struct();
Fig10Data.N = N;
Fig10Data.K = K;
Fig10Data.p = p;
Fig10Data.networkSeed = networkSeed;
Fig10Data.numRepeat = numRepeat;
Fig10Data.h = params0.h;
Fig10Data.T = params0.T;
Fig10Data.epsH = epsH;
Fig10Data.sigmaPresent = sigmaPresent;
Fig10Data.sigmaRef54 = sigmaRef54;
Fig10Data.GcPresent = GcPresent;
Fig10Data.GcRef54 = GcRef54;
Fig10Data.Present = Present;
Fig10Data.Ref54 = Ref54;

save('Fig10.mat','Fig10Data','-v7.3')
toc


end
function S = runMethod(mode,Gc,x0Cell,params0)

numG = numel(Gc);
numRepeat = numel(x0Cell);
M = params0.M_edge;
N = params0.N;
Results = NaN(numG*numRepeat,3);

parfor job = 1:numG*numRepeat
    [gIdx,rIdx] = ind2sub([numG,numRepeat],job);

    params = params0;
    params.mode = mode;
    params.gc = Gc(gIdx);

    if strcmp(mode,'present')
        params.q0 = ones(M,1);
    else
        params.k0 = Gc(gIdx)*ones(N,1);
    end

    R = FHNres_Ref54Present(x0Cell{rIdx},params);

    Results(job,:) = [R.hSF,R.HJ,R.meanCoupling];
end

All = reshape(Results,numG,numRepeat,3);

S.HSF.mean = squeeze(mean(All(:,:,1),2,'omitnan'));
S.HSF.std = squeeze(std(All(:,:,1),0,2,'omitnan'));
S.HJ.mean = squeeze(mean(All(:,:,2),2,'omitnan'));
S.HJ.std = squeeze(std(All(:,:,2),0,2,'omitnan'));
S.MeanCoupling.mean = squeeze(mean(All(:,:,3),2,'omitnan'));
S.MeanCoupling.std = squeeze(std(All(:,:,3),0,2,'omitnan'));
S.all = All;

end
