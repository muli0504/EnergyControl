function FigA2
clc
close all
tic
addpath('NetworkAnalyse')
beta = (0.01:0.001:1.0).';
LE = NaN(numel(beta),3);
if isempty(gcp('nocreate'))
    parpool;
end
parfor i = 1:numel(beta)
    LE(i,:) = LyapunovBeta(beta(i));
end
FigA2Data = struct();
FigA2Data.beta = beta;
FigA2Data.LE = LE;
FigA2Data.LLE = LE(:,1);
save('FigA2.mat','FigA2Data','-v7.3')
toc

end
