% sweep parameters of SPP model and plot distribution of peak delay times
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(2:0.5:7);
betaValues = 2.^(0:0.5:7);
numRepeats = 10;
T = 1000;
burnIn = 500;
N = 100;
L = 2;
r0 = 1;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

binWidth = 3;
maxLag = 20;
lagValues = -27:27; % should allow for an integer numbe of bins
% use two more bins than wanted (one for nicer plotting, one to discard edges)
nLagValues = length(lagValues);
minCorr = 0.5;
minOrder = 0.1; % as computing correlations for all pairs and delays is expensive, only do for relatively ordered collectives

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',12,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

precision = 2;

%% load results
binCentres = (min(lagValues)+2*binWidth):binWidth:(max(lagValues)-2*binWidth);
histPeakLag = NaN(numAlphas,numBetas,numRepeats,(nLagValues - 1)/binWidth - 1);
histPeakMeanAbs = NaN(numAlphas,numBetas,numRepeats);
histPeakVar = NaN(numAlphas,numBetas,numRepeats);
histPeakKurt = NaN(numAlphas,numBetas,numRepeats);
histPeakKurt0 = NaN(numAlphas,numBetas,numRepeats);
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        distributionFig = figure;
        distributionFig.Color='none';
        hold on
        plotCurrent = 0;
        order = NaN(numRepeats,1);
        for repCtr = 1:numRepeats
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ... %'_selfAlign' ...
                '_run' num2str(repCtr) '.mat'];
            out = load(filename);
            % discard burn-in
            out.cells = out.cells(:,:,burnIn:end);
            order(repCtr) = mean(orderParameter(out.cells));
            if order(repCtr) >= minOrder
                dirCrossCorr = NaN(N*(N-1)/2,nLagValues);
                % only need to go over each pair once, since Cij(tau) = Cji(-tau)
                
                % calculate cross-correlations
                for lagCtr = 1:nLagValues % can be parfor-ed
                    lag = lagValues(lagCtr);
                    dirCrossCorr(:,lagCtr) = mean(directionalCrossCorrelation(out.cells,lag,r0),2);
                end
                
                % only keep cross-correlations if cells have been neighbours (within r0)
                pairDistances = NaN(N*(N-1)/2,T-burnIn);
                for timeCtr=1:T-burnIn % can be parfor-ed
                    pairDistances(:,timeCtr) = pdist(out.cells(:,1:3,timeCtr));
                end
                notNeighbours = min(pairDistances,[],2)>r0;
                dirCrossCorr(notNeighbours,:) = NaN;
                
                % find and plot peaks in crosscorr curves
                peakLags = NaN(N,N);
                peakCorrs = NaN(N,N);
                
                ij = 0;
                for jj = 1:(N-1)
                    for ii = (jj+1):N
                        ij = ij+1;
                        % if the two cells were never neighbours, their crossCorr = NaN
                        if any(~isnan(dirCrossCorr(ij,:)))
                            [pks, locs] = findpeaks(dirCrossCorr(ij,:),lagValues,...
                                'NPeaks',1,'SortStr','descend');
                            if ~isempty(pks), peakCorrs(ii,jj) = pks; end
                            if ~isempty(locs), peakLags(ii,jj) = locs; end
                            peakCorrs(jj,ii) = peakCorrs(ii,jj); % correlation is the same for the symmetric pair
                            peakLags(jj,ii) = -peakLags(ii,jj);% Cij(tau) = Cji(-tau)
                        end
                    end
                end
                
                % only plot lag times for corr>=minCorr - don't use
                % 'MinPeakHeight' parameter in findpeaks as it throws
                % warnings
                threshPeaks = peakCorrs>=minCorr;
                % discard peak lag times outside relevant range from
                % further calculation
                threshPeaks(abs(peakLags(threshPeaks))>maxLag) = false;
                if any(threshPeaks(:))
                    plotCurrent = 1;
                    %% save distribution
                    histPeakLag(alphaCtr,betaCtr,repCtr,:) = histcounts(peakLags(threshPeaks),...
                        'BinWidth',binWidth,'BinLimits',[min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                        'Normalization','probability');
                    histPeakMeanAbs(alphaCtr,betaCtr,repCtr) = mean(abs(peakLags(threshPeaks)));
                    histPeakVar(alphaCtr,betaCtr,repCtr) = var(peakLags(threshPeaks));
                    histPeakKurt(alphaCtr,betaCtr,repCtr) = kurtosis(peakLags(threshPeaks),1);
                    histPeakKurt0(alphaCtr,betaCtr,repCtr) = kurtosis(peakLags(threshPeaks),0);
                end
            end
        end
        if plotCurrent
            boundedline(binCentres,...
                squeeze(nanmean(histPeakLag(alphaCtr,betaCtr,:,2:end-1),3)),...
                squeeze(nanstd(histPeakLag(alphaCtr,betaCtr,:,2:end-1),0,3)),'.-');
            ylim([0 max(max(histPeakLag(alphaCtr,betaCtr,:,:)))])
            xlim([-maxLag maxLag])
            % format figure
            ax = gca;
            ax.Box = 'off';
            ax.XLabel.String = '\tau_C';
            ax.YLabel.String = 'P';
            ax.Title.String = ['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3),...
                ', \langle\Phi\rangle = ' num2str(nanmean(order),2)];
            ax.Title.FontWeight = 'normal';
            %% export figure
            filename = ['manuscript/figures/diagnostics/delayDist_T' num2str(T) '_N' num2str(N) ...
                '_L' num2str(L) '_a' num2str(alpha,precision) ...
                '_b' num2str(beta,precision) ... %'_selfAlign'
                ];
            %export_fig([filename '.pdf'],'-depsc')
            set(distributionFig,'PaperUnits','centimeters')
            exportfig(distributionFig,[filename '.eps'],exportOptions);
            system(['epstopdf ' filename '.eps']);
            system(['rm ' filename '.eps']);
        end
        close(distributionFig)
    end
end
%% plot multi-line diagram of variance and kurtosis
exportOptions.FontSize = 14;
exportOptions.Width = '15';
legendString = num2str(round(alphaValues(1:2:end))');
absFig = figure;
boundedline(betaValues,nanmean(histPeakMeanAbs(1:2:end,:,:),3),...
    permute(nanstd(histPeakMeanAbs(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakMeanAbs(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
ax1 = gca;
ax1.XScale = 'log';
ax1.XLabel.String = '\beta'; ax1.YLabel.String = '\langle|\tau_C|\rangle';
absFig.Color='none'; ax1.Box = 'on';
xlim([1 max(betaValues)])
ylim([0 max(max(max(histPeakMeanAbs)))])
% make an inset of non-log scale
inset = axes('position',[0.25 0.5 0.4 0.4]);
boundedline(betaValues,nanmean(histPeakMeanAbs(1:2:end,:,:),3),...
    permute(nanstd(histPeakMeanAbs(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakMeanAbs(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
inset.XLim = [1 max(betaValues)];
inset.YLim = [0 max(max(histPeakMeanAbs(1:2:end,end,:,:)))];
inset.Box = 'on';
legH = legend(ax1,legendString,'Location','NorthEast');
legH.Title.String = '\alpha';
% save figure
filename = ['manuscript/figures/meanabsDelayDiagram_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L) ... %'_selfAlign'...
    ];
set(absFig,'PaperUnits','centimeters')
exportfig(absFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']);

sigmaFig = figure;
boundedline(betaValues,nanmean(sqrt(histPeakVar(1:2:end,:,:)),3),...
    permute(nanstd(sqrt(histPeakVar(1:2:end,:,:)),0,3)./...
    sqrt(sum(~isnan(histPeakVar(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
ax1 = gca;
ax1.XScale = 'log';
ax1.XLabel.String = 'attraction-repulsion strength \beta'; ax1.YLabel.String = 'heterogeneity \sigma(\tau_C)';
absFig.Color='none'; ax1.Box = 'on';
xlim([1 max(betaValues)])
ylim([0 1+ceil(max(max(max(sqrt(histPeakVar(1:2:end,:,:))))))])
% make an inset of non-log scale
inset = axes('position',...
    [0.33 0.59 0.4 0.335]...
    ...%[0.43 0.1275 0.385 0.265] % selfalign
);
boundedline(betaValues,nanmean(sqrt(histPeakVar(1:2:end,:,:)),3),...
    permute(nanstd(sqrt(histPeakVar(1:2:end,:,:)),0,3)./...
    sqrt(sum(~isnan(histPeakVar(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
inset.XLim = [1 max(betaValues)];
inset.YLim = [4 12.5];
inset.YTick = [4 8 12];
inset.Box = 'on';
%inset.XAxisLocation = 'top'; %selfalign
% legH = legend(ax1,legendString,'Location','SouthEast');
% legH.Title.String = '\alpha';
% legH.Position = [0.7325 0.17 0.1 0.199]...%[0.805 0.185 0.12 0.199] %selfalign
% ;
% save figure
filename = ['manuscript/figures/varianceDelayDiagram_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L) ... %'_selfAlign'...
    ];
set(sigmaFig,'PaperUnits','centimeters')
exportfig(sigmaFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

kurtFig = figure;
boundedline(betaValues,nanmean(histPeakKurt(1:2:end,:,:),3),...
    permute(nanstd(histPeakKurt(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakKurt(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
ax1 = gca;
ax1.XScale = 'log';
ax1.XLabel.String = '\beta'; ax1.YLabel.String = 'Kurt(\tau_C)';
absFig.Color='none'; ax1.Box = 'on';
xlim([1 max(betaValues)])
ylim([0 max(max(max(histPeakKurt(1:2:end,:,:))))])
% make an inset of non-log scale
inset = axes('position',[0.25 0.5 0.4 0.4]);
boundedline(betaValues,nanmean(histPeakKurt(1:2:end,:,:),3),...
    permute(nanstd(histPeakKurt(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakKurt(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
inset.XLim = [1 max(betaValues)];
inset.YLim = [0 max(max(histPeakKurt(1:2:end,end,:,:)))];
inset.Box = 'on';
legH = legend(ax1,legendString,'Location','NorthEast');
legH.Title.String = '\alpha';
% save figure
filename = ['manuscript/figures/kurtosisDelayDiagram_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L) ... %'_selfAlign'...
    ];
set(kurtFig,'PaperUnits','centimeters')
exportfig(kurtFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

kurt0Fig = figure;
boundedline(betaValues,nanmean(histPeakKurt0(1:2:end,:,:),3),...
    permute(nanstd(histPeakKurt0(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakKurt0(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
ax1 = gca;
ax1.XScale = 'log';
ax1.XLabel.String = '\beta'; ax1.YLabel.String = 'Kurt_0(\tau_C)';
absFig.Color='none'; ax1.Box = 'on';
xlim([1 max(betaValues)])
ylim([0 max(max(max(histPeakKurt0(1:2:end,:,:))))])
% make an inset of non-log scale
inset = axes('position',[0.25 0.5 0.4 0.4]);
boundedline(betaValues,nanmean(histPeakKurt0(1:2:end,:,:),3),...
    permute(nanstd(histPeakKurt0(1:2:end,:,:),0,3)./...
    sqrt(sum(~isnan(histPeakKurt0(1:2:end,:,:)),3)),[3 2 1]),'alpha','.-','nan','gap')
inset.XLim = [1 max(betaValues)];
inset.YLim = [0 max(max(histPeakKurt0(1:2:end,end,:,:)))];
inset.Box = 'on';
legH = legend(ax1,legendString,'Location','NorthEast');
legH.Title.String = '\alpha';
% save figure
filename = ['manuscript/figures/kurtosisDebiasedDelayDiagram_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L) ... %'_selfAlign'...
    ];
set(kurt0Fig,'PaperUnits','centimeters')
exportfig(kurt0Fig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);
%% surface plot
figure
surface(betaValues,alphaValues,mean(sqrt(histPeakVar),3),...
    'FaceColor','interp','EdgeColor','none')
xlabel('\beta'),ylabel('\alpha')
set(gca,'xscale','log','yscale','log')
xlim([1 128]),ylim([4 128])
maxSig = round(max(max(mean(sqrt(histPeakVar),3))));
colormap(flipud(gray(length(0:0.5:maxSig) - 1)))
cb = colorbar;
cb.Label.String = '\sigma(\tau_c)';
cb.Label.Rotation = 0;
cb.Label.HorizontalAlignment = 'left';
caxis([0 maxSig])
box on
% save figure
filename = ['manuscript/figures/varianceDelayDiagramSurface_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L) ... %'_selfAlign'...
    ];
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);
%%
save(['delayCorrResults' ... %'_selfAlign'
    '.mat'],'alphaValues','betaValues','histPeakMeanAbs',...
    'histPeakVar','histPeakKurt0','histPeakKurt')