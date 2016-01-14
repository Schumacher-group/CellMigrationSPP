% sweep parameters of SPP model to generate simulation data
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

alphaValues = {2.^(1:7); 4:8; 4:8};
betaValues = {2.^(-1:7); 2.^(1:0.5:6); 2.^(1:0.5:6)};
numRepeats = 10;
T = 1000;
Nvalues = [100, 20, 10];
Nlvalues = [10, 2, 1];
Lvalues = [2, 1.1, 0.9];

bcs = 'free';
precision = 2;
rng('shuffle')
for permCtr = 1:length(Nvalues)
    N = Nvalues(permCtr);
    Nl = Nlvalues(permCtr);
    L = Lvalues(permCtr);
    numAlphas = length(alphaValues{permCtr});
    numBetas = length(betaValues{permCtr});
    for selfAlign = [false] %, true]
        for repCtr = randsample(1:numRepeats,numRepeats,0)
            for alphaCtr = randsample(1:numAlphas,numAlphas,0)
                alpha = alphaValues{permCtr}(alphaCtr);
                for betaCtr = randsample(1:numBetas,numBetas,0)
                    beta = betaValues{permCtr}(betaCtr);
                    filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                        '_Nl' num2str(Nl) '_L' num2str(L,precision) '_' bcs ...
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