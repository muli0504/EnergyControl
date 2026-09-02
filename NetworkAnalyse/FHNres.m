function R = FHNres(x0,params)
h = params.h;
N = params.N;
T = params.T;
M = params.M_edge;
Amp = params.Amp(:);
omega = 2*pi*params.freq(:);
transient = params.transient;
K = T-transient;
calcCSF = getOpt(params,'calcCSF',false);
calcHSF = getOpt(params,'calcHSF',false);
calcEn = getOpt(params,'calcEn',false);
calcECF = getOpt(params,'calcECF',false);
calcHJ = getOpt(params,'calcHJ',false);
calcGStat = getOpt(params,'calcGStat',false);
calcSync = getOpt(params,'calcEnergySyncTime',false);
calcMap = getOpt(params,'calcFig7Map',false);
trackPairs = getOpt(params,'edgeTrackPairs',[]);
trackStep = max(1,round(getOpt(params,'edgeTrackStep',10)));
mapStep = max(1,round(getOpt(params,'fig7MapStep',10)));
threshold = getOpt(params,'energySyncThreshold',1e-3);
holdWindow = max(1,round(getOpt(params,'energySyncHoldWindow',1000)));
startStep = max(1,round(getOpt(params,'energySyncStartStep',1)));
if ~isequal(size(x0),[N,3])
    error('FHNres: x0 must be N-by-3.')
end
q0 = getOpt(params,'q0',ones(M,1));
q0 = q0(:);
if numel(q0) ~= M
    error('FHNres: q0 length mismatch.')
end
Y = [x0(:);q0];
csf = initStat(N);
hsf = initStat(N);
En = zeros(N,1);
ecf1 = 0;
ecf2 = 0;
HJsum = 0;
gMeanSum = 0;
gCVsum = 0;
syncTime = NaN;
syncReached = false;
syncCount = 0;
syncCandidate = NaN;
syncErrMin = Inf;
syncErrFinal = NaN;
if ~isempty(trackPairs)
    [trackIdx,trackFound] = findPairs(params.edge_i,params.edge_j,trackPairs);
    nTrack = floor((K-1)/trackStep)+1;
    trackTime = NaN(nTrack,1);
    trackG = NaN(nTrack,size(trackPairs,1));
    trackCounter = 1;
else
    trackIdx = [];
    trackFound = false(0,1);
    trackTime = [];
    trackG = [];
    trackCounter = 1;
end
if calcMap
    nMap = floor((K-1)/mapStep)+1;
    mapTime = NaN(nMap,1);
    VmMap = NaN(nMap,N,'single');
    EnergyMap = NaN(nMap,N,'single');
    edgeHJMap = NaN(nMap,M,'single');
    edgeGMap = NaN(nMap,M,'single');
    mapCounter = 1;
else
    mapTime = [];
    VmMap = [];
    EnergyMap = [];
    edgeHJMap = [];
    edgeGMap = [];
    mapCounter = 1;
end
t = 0;
for step = 1:T
    I1 = Amp.*cos(omega*t);
    I2 = Amp.*cos(omega*(t+0.5*h));
    I3 = I2;
    I4 = Amp.*cos(omega*(t+h));
    K1 = FHNfunction(Y,I1,params);
    K2 = FHNfunction(Y+0.5*h*K1,I2,params);
    K3 = FHNfunction(Y+0.5*h*K2,I3,params);
    K4 = FHNfunction(Y+h*K3,I4,params);
    Y = Y + h*(K1+2*K2+2*K3+K4)/6;
    q = protectQ(Y(3*N+1:3*N+M),params);
    Y(3*N+1:3*N+M) = q;
    X = reshape(Y(1:3*N),N,3);
    t = t+h;
    if calcSync && step >= startStep
        H = NodeEnergy(X,I4,params);
        err = sqrt(mean((H-mean(H,'omitnan')).^2,'omitnan'));
        syncErrFinal = err;
        if isfinite(err)
            syncErrMin = min(syncErrMin,err);
            if err <= threshold
                if syncCount == 0
                    syncCandidate = t;
                end
                syncCount = syncCount+1;
                if ~syncReached && syncCount >= holdWindow
                    syncTime = syncCandidate;
                    syncReached = true;
                end
            else
                syncCount = 0;
                syncCandidate = NaN;
            end
        else
            syncCount = 0;
            syncCandidate = NaN;
        end
    end
    if step <= transient
        continue
    end
    localStep = step-transient;
    if ~isempty(trackPairs) && mod(localStep-1,trackStep)==0
        g = Gedge_from_q(q,params);
        trackTime(trackCounter) = t;
        ok = trackFound(:);
        trackG(trackCounter,ok) = g(trackIdx(ok));
        trackCounter = trackCounter+1;
    end
    if calcCSF
        csf = updateStat(csf,X(:,1),K);
    end
    if calcHJ || calcGStat
        g = Gedge_from_q(q,params);
        if calcHJ
            dx = X(params.edge_i,1)-X(params.edge_j,1);
            HJsum = HJsum + mean(g.*dx.^2);
        end
        if calcGStat
            gMean = mean(g);
            gStd = std(g,0);
            gMeanSum = gMeanSum + gMean;
            gCVsum = gCVsum + gStd/(gMean+eps);
        end
    end
    needH = calcHSF || calcEn || calcECF || calcMap;
    if needH
        H = NodeEnergy(X,I4,params);
        if calcHSF
            hsf = updateStat(hsf,H,K);
        end
        if calcEn
            En = En + H/K;
        end
        if calcECF
            rE = mean(H)^2/(mean(H.^2)+eps);
            ecf1 = ecf1 + rE^2/K;
            ecf2 = ecf2 + rE/K;
        end
        if calcMap && mod(localStep-1,mapStep)==0
            g = Gedge_from_q(q,params);
            dx = X(params.edge_i,1)-X(params.edge_j,1);
            mapTime(mapCounter) = t;
            VmMap(mapCounter,:) = single(X(:,1)).';
            EnergyMap(mapCounter,:) = single(H).';
            edgeHJMap(mapCounter,:) = single(g.*dx.^2).';
            edgeGMap(mapCounter,:) = single(g).';
            mapCounter = mapCounter+1;
        end
    end
end
R.CSF = NaN;
R.hSF = NaN;
R.En = NaN(N,1);
R.ecf = NaN;
R.HJ = NaN;
R.gEdge_mean_time = NaN;
R.gEdge_cv_time = NaN;
if calcCSF
    R.CSF = syncFactor(csf);
end
if calcHSF
    R.hSF = syncFactor(hsf);
end
if calcEn
    R.En = En;
end
if calcECF
    R.ecf = N*max(ecf1-ecf2^2,0);
end
if calcHJ
    R.HJ = HJsum/K;
end
if calcGStat
    R.gEdge_mean_time = gMeanSum/K;
    R.gEdge_cv_time = gCVsum/K;
end
R.energySyncTime = syncTime;
R.energySyncReached = syncReached;
R.energySyncErrMin = syncErrMin;
R.energySyncErrFinal = syncErrFinal;
if ~isempty(trackPairs)
    last = trackCounter-1;
    R.edgeTrackPairs = trackPairs;
    R.edgeTrackFound = trackFound;
    R.edgeTrackTime = trackTime(1:last);
    R.edgeTrackG = trackG(1:last,:);
else
    R.edgeTrackPairs = zeros(0,2);
    R.edgeTrackFound = false(0,1);
    R.edgeTrackTime = [];
    R.edgeTrackG = [];
end
R.edgePairs = [params.edge_i(:),params.edge_j(:)];
if calcMap
    last = mapCounter-1;
    R.mapTime = mapTime(1:last);
    R.VmMap = VmMap(1:last,:);
    R.EnergyMap = EnergyMap(1:last,:);
    R.edgeHJMap = edgeHJMap(1:last,:);
    R.edgeGMap = edgeGMap(1:last,:);
else
    R.mapTime = [];
    R.VmMap = [];
    R.EnergyMap = [];
    R.edgeHJMap = [];
    R.edgeGMap = [];
end
end

function v = getOpt(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = defaultValue;
end
end

function q = protectQ(q,params)
q = max(q(:),getOpt(params,'qMin',1e-12));
if isfield(params,'qMax') && ~isempty(params.qMax)
    q = min(q,params.qMax);
end
if getOpt(params,'rescaleQMean',false)
    qMean = mean(q);
    if ~isfinite(qMean) || qMean <= 0
        error('FHNres: invalid q mean.')
    end
    q = q/qMean;
end
end

function S = initStat(N)
S.s1 = 0;
S.s2 = 0;
S.s3 = zeros(N,1);
S.s4 = zeros(N,1);
end

function S = updateStat(S,x,K)
F = mean(x);
S.s1 = S.s1 + F^2/K;
S.s2 = S.s2 + F/K;
S.s3 = S.s3 + x.^2/K;
S.s4 = S.s4 + x/K;
end

function sf = syncFactor(S)
den = mean(S.s3-S.s4.^2);
if den > 0
    sf = (S.s1-S.s2^2)/den;
else
    sf = NaN;
end
end

function [idx,found] = findPairs(i,j,pairs)
pairs = double(pairs);
idx = NaN(size(pairs,1),1);
found = false(size(pairs,1),1);
for k = 1:size(pairs,1)
    a = min(pairs(k,:));
    b = max(pairs(k,:));
    n = find(i==a & j==b,1);
    if ~isempty(n)
        idx(k) = n;
        found(k) = true;
    end
end
end
