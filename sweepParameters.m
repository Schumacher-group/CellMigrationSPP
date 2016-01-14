% sweep parameters of SPP model to calculate order phase diagrams
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

alphaValues = 2.^(2:3);
betaValues = 2.^(0:4);
T = 600;
burnIn = 500;
N = 100;
L = 2;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

bcs = {'noflux','noflux','free'};

precision = 3;
trajectoryFig = figure;

for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        % run simulation
        cells = runSPP(T,N,L,alpha,beta,'bc',bcs);
        % calculate order
        order = orderParameter(cells);
        %plot trajectories
        %     set(0,'CurrentFigure',trajectoryFig)
        subplot(numAlphas,numBetas,(alphaCtr - 1)*numBetas + betaCtr)
        plot3(squeeze(cells(:,1,burnIn:end))',squeeze(cells(:,2,burnIn:end))',...
            squeeze(cells(:,3,burnIn:end))','-','LineWidth',1)
        axis image
        ax = gca;
        ax.Title.String = ['\alpha = ' num2str(alpha,precision) ', \beta = ' num2str(beta,precision)...
            ', \Phi = ' num2str(mean(order(burnIn:end)),precision)];
        ax.Title.FontWeight = 'normal';
        ax.XLabel.String = 'x';
        ax.YLabel.String = 'y';
        ax.ZLabel.String = 'z';
        %     set(0,'CurrentFigure',dzFig)
        %     subplot(numAlphas,numBetas,(alphaCtr - 1)*numBetas + betaCtr)
        %     plot(squeeze(cells(:,3,burnIn:end))')
        %     % plot order
        %     set(0,'currentfigure',orderFig)
        %     plot(order)
        %     xlabel('time'), ylabel('relative order parameter')
    end
end