% sweep parameters of SPP model to calculate order phase diagrams
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

alphaValues = 2.^(1:3);
betaValues = 2.^(0:4);
numRepeats = 10;
T = 1000;
N = 9;
L = 0.6;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

bcs = {'free','noflux','noflux'};
precision = 2;
selfAlign = false;
rng('shuffle')
for repCtr = randsample(1:numRepeats,numRepeats,0)
    for alphaCtr = randsample(1:numAlphas,numAlphas,0)
        alpha = alphaValues(alphaCtr);
        for betaCtr = randsample(1:numBetas,numBetas,0)
            beta = betaValues(betaCtr);
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) '_' bcs{1} '-' bcs{2} '-' bcs{3} ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...%'_selfAlign' ...
                '_run' num2str(repCtr) '.mat'];
            if ~exist(filename,'file')
                % run simulation
                rng(repCtr) % each run from a different parameter combination should have the same random numbers
                cells = runSPP(T,N,L,alpha,beta,'bc',bcs,'selfAlign',selfAlign);
                % only keep position, not directional data
                cells = cells(:,[x y z],:);
                save(filename,'cells','T','N','L','alpha','beta','selfAlign','bcs')
            end
        end
    end
end