% sweep parameters of SPP model to generate simulation data
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

alphaValues = 2.^[2:3];
betaValues = 2.^[0:4];
numRepeats = 10;
T = 1000;
Nvalues = [100 30 10];
Lvalues = [2 sqrt(3/10)*2 sqrt(1/10)*2];
Nlvalues = [10, 3, 1];
numAlphas = length(alphaValues);
numBetas = length(betaValues);

bcs = {'free', 'noflux', 'noflux'};
precision = 2;
rng('shuffle')
for permCtr = 1:length(Nvalues)
    N = Nvalues(permCtr);
    Nl = Nlvalues(permCtr);
    L = Lvalues(permCtr);
    for selfAlign = [false] %, true]
        for repCtr = randsample(1:numRepeats,numRepeats,0)
            for alphaCtr = randsample(1:numAlphas,numAlphas,0)
                alpha = alphaValues(alphaCtr);
                for betaCtr = randsample(1:numBetas,numBetas,0)
                    beta = betaValues(betaCtr);
                    filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                        '_Nl' num2str(Nl) '_L' num2str(L,precision) '_' bcs{1} '-' bcs{2} '-' bcs{3} ...
                        '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
                        '_selfAlign' num2str(selfAlign) '_run' num2str(repCtr) '.mat'];
                    if ~exist(filename,'file')
                        % run simulation
                        rng(repCtr) % each run from a different parameter combination should have the same random numbers
                        cells = runSPP(T,N,L,alpha,beta,'bc',bcs,...
                            'Nl',Nl,'selfAlign',selfAlign);
                        % only keep position, not directional data
                        cells = cells(:,[x y z],:);
                        save(filename,'cells','T','N','Nl','L','alpha','beta','selfAlign','bcs')
                    end
                end
            end
        end
    end
end