 % compare the distributions of peak delay times for different parameters
 % under confinement
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^[2:3];
betaValues = 2.^[0:4];
numRepeats = 10;
T = 1000;
burnIn = 500;
NValues = [100 25 9];
LValues = [2 1 0.6];
% for plotting results from the model with leaders uncomment below
% NValues = [100 30 10];
% LValues = [2 sqrt(3/10)*2 sqrt(1/10)*2];
% Nlvalues = [10, 3, 1];
r0 = 1;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

plotColors = lines(numBetas*numAlphas);
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',9,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

precision = 2;
%% load reference simulation results (no confinement)
reference = load('delayCorrResults.mat','alphaValues','betaValues','histPeakVar');
% for plotting results from the model with leaders uncomment below
% reference = load('delayCorrResults_Nl10_selfAlign0.mat','alphaValues','betaValues','histPeakVar');

%% load results to compare with
for permCtr = 1:3
    L = LValues(permCtr);
    N = NValues(permCtr);
       results(permCtr) = load(['delayCorrResultsConfined_N' num2str(N) '_L' num2str(L,precision) '.mat'],...        
'alphaValues','betaValues','histPeakVar');
    % for plotting results from the model with leaders uncomment below
%     Nl = Nlvalues(permCtr);
%     results(permCtr) = load(['delayCorrResultsConfined_N' num2str(N) '_Nl' num2str(Nl) '_L' num2str(L,precision) '.mat'],...
end
%% plot distribution widths
numPlots = 3*length(alphaValues)*length(betaValues);
plotColors = lines(numPlots);
plotCtr = 1;
figure, hold on
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        plotSeries = sqrt(squeeze(reference.histPeakVar(reference.alphaValues==alpha,reference.betaValues==beta,:)));
        for permCtr = 1:3
            plotSeries = [plotSeries, sqrt(squeeze(results(permCtr).histPeakVar(alphaCtr,betaCtr,:)))];
        end

        plot(nanmean(plotSeries),'.-','MarkerSize',12,'LineWidth',2,'Color',...
        [plotColors(plotCtr,:) 0.25],'MarkerFaceColor',plotColors(plotCtr,:))
        plotCtr = plotCtr + 1;
    end
end
%% annotate plot
ax = gca;
ax.Box = 'on';
ax.XTick = 1:4;
ax.XTickLabel = {'free','N=100, L=2','N=25, L=1','N=9, L=0.6'};
    % for plotting results from the model with leaders uncomment below
% ax.XTickLabel = {'free','N=100, L=2','N=30, L=1.1','N=10, L=0.63'};
ax.XLabel.String = 'confinement';
ax.YLabel.String = 'heterogeneity \sigma(\tau_c)';
legH.Title.String = '\beta';


%% export figure
set(gcf,'PaperUnits','centimeters')
filename = ['manuscript/figures/delayVarConfinement'];
    % for plotting results from the model with leaders uncomment below
% filename = ['manuscript/figures/delayVarConfinementWithLeaders'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);
