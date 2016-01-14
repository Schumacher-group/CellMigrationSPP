% sweep parameters of SPP model to calculate order phase diagrams
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

alphaValues = 2.^(-1:4);
betaValues = 2.^(-2:2:8);
numRepeats = 1;
T = 1000;
N = 10;
L = [1, 1, 0.05];
numAlphas = length(alphaValues);
numBetas = length(betaValues);

bcs = {'free', 'free', 'noflux'};

precision = 2;
selfAlign = false;

for repCtr = 1:numRepeats;
    for alphaCtr = 1:numAlphas
        alpha = alphaValues(alphaCtr);
        for betaCtr = 1:numBetas
            beta = betaValues(betaCtr);
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L(1),precision) '_' num2str(L(2),precision)...
                '_' num2str(L(3),precision) '_' bcs{1} '-' bcs{2} '-' bcs{3}...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision)...
                '_run' num2str(repCtr) '.mat'];
% %             if ~exist(filename,'file')
                % run simulation
                rng(repCtr) % each run from a different parameter combination should have the same random numbers
                cells = runSPP(T,N,L,alpha,beta,'bc',bcs);
                % only keep position, not directional data
                cells = cells(:,[x y z],:);
                save(filename,'cells','T','N','L','alpha','beta','selfAlign','bcs')
% %             end
        end
    end
end