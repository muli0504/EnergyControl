function Fig2
clc
close all
tic
addpath('NetworkAnalyse')
fList = (0:0.0001:0.25).';
if isempty(gcp('nocreate'))
    parpool;
end
LLE = NaN(numel(fList),1);
MeanEnergy = NaN(numel(fList),1);
parfor i = 1:numel(fList)
    ly = LyapunovF(fList(i));
    LLE(i) = ly(1);
    MeanEnergy(i) = SingleEnergyF(fList(i));
end
rng(1,'twister')
initXYZ = 2*rand(200,3);
[LEs,InitialMeanEnergy,boundedFlag] = InitialConditionLLE(initXYZ);
Fig2Data = struct();
Fig2Data.f = fList;
Fig2Data.LLE = LLE;
Fig2Data.MeanEnergy = MeanEnergy;
Fig2Data.initXYZ = initXYZ;
Fig2Data.InitialLE = LEs;
Fig2Data.InitialLLE = LEs(:,1);
Fig2Data.InitialMeanEnergy = InitialMeanEnergy;
Fig2Data.boundedFlag = boundedFlag;
Fig2Data.chaosFlag = boundedFlag & isfinite(LEs(:,1)) & LEs(:,1)>0;
save('Fig2.mat','Fig2Data','-v7.3')
toc

end
