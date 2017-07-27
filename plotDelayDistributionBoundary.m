% plot distribution of peak delay times for parts of cell population at the
% boundary and centre
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 4;
betaValues = 4;
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
% use two more bins than wanted (one fpr nicer plotting, one to discard edges)
nLagValues = length(lagValues);
minCorr = 0.5;
minOrder = 0.1; % as computing correlations for all pairs and delays is expensive, only do for relatively ordered collectives

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

precision = 2;

%% load results
binCentres = (min(lagValues)+2*binWidth):binWidth:(max(lagValues)-2*binWidth);

for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        distributionFig = figure;
        distributionFig.Color='none';
        hold on
        histPeakLagBoundary = NaN(numRepeats,(nLagValues - 1)/binWidth - 1);
        histPeakLagCore = NaN(numRepeats,(nLagValues - 1)/binWidth - 1);
        histPeakLagMixed = NaN(numRepeats,(nLagValues - 1)/binWidth - 1);
        for repCtr = 1:numRepeats
            plotCurrent = 0;
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
                '_run' num2str(repCtr) '.mat'];
            out = load(filename);
            % discard burn-in
            out.cells = out.cells(:,:,burnIn:end);
            order = mean(orderParameter(out.cells));
            if order >= minOrder
                dirCrossCorr = NaN(N*(N-1)/2,nLagValues);
                % only need to go over each pair once, since Cij(tau) = Cji(-tau)
                
                % calculate cross-correlations
                for lagCtr = 1:nLagValues % can be parfor-ed
                    lag = lagValues(lagCtr);
                    dirCrossCorr(:,lagCtr) = mean(directionalCrossCorrelation(out.cells,lag,r0),2);
                end
                
                % only keep cross-correlations if cells have been neighbours (within r0)
                pairDistances = NaN(N*(N-1)/2,T-burnIn);
                boundaryMtx = false(N,T-burnIn);
                for timeCtr=1:T-burnIn % can be parfor-ed
                    pairDistances(:,timeCtr) = pdist(out.cells(:,1:3,timeCtr));
                    % check which cells are at the boundary
                    [~, boundaryPointIDs] = getBoundaryPoints(out.cells(:,1:3,timeCtr),r0);
                    boundaryMtx(boundaryPointIDs,timeCtr) = true;
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
                    % plot separate distributions for cells at boundary,
                    % core, and neither
                    boundaryCells = mean(boundaryMtx,2)>0.9; % define boundary cells as those that are at the boundary more than 90% of the time
                    coreCells = mean(boundaryMtx,2)<0.1; % define core cells as those at the boundary less than 10% of the time
                    histPeakLagBoundary(repCtr,:) = histcounts(peakLags(threshPeaks(boundaryCells,:)),...
                        'BinWidth',binWidth,'BinLimits',...
                        [min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                        'Normalization','probability');
                    histPeakLagCore(repCtr,:) = histcounts(peakLags(threshPeaks(coreCells,:)),...
                        'BinWidth',binWidth,'BinLimits',...
                        [min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                        'Normalization','probability');
                    histPeakLagMixed(repCtr,:) = histcounts(peakLags(threshPeaks(~(coreCells|boundaryCells),:)),...
                        'BinWidth',binWidth,'BinLimits',...
                        [min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                        'Normalization','probability');
                end
            end
        end
        if plotCurrent
            boundedline(binCentres,[nanmean(histPeakLagMixed(:,2:end-1));...
                nanmean(histPeakLagCore(:,2:end-1)); ...
                nanmean(histPeakLagBoundary(:,2:end-1))],... % averages
                permute(repmat([nanstd(histPeakLagMixed(:,2:end-1));...
                nanstd(histPeakLagCore(:,2:end-1));...
                nanstd(histPeakLagBoundary(:,2:end-1))],1,1,2),[2 3 1]),... % stdev
                'alpha','.-','cmap',[0.5*[1 1 1]; 0.7 0 1; 1 0.7 0]);
            % plotting multiple bounded lines in one go works better
            legend('mixed','core','boundary')
            % format figure
            ax = gca;
            ax.XLim = [-maxLag maxLag];
            ax.YLim = [0 max(max([histPeakLagBoundary;
                histPeakLagCore; histPeakLagMixed]))];
            ax.Box = 'on';
            ax.XLabel.String = 'peak delay time \tau_C';
            ax.YLabel.String = 'probability P';
            ax.YTickLabel = num2str(ax.YTick','%1.2f');
            ax.Title.String = ['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3)];
            ax.Title.FontWeight = 'normal';
            %% export figure
            filename = ['manuscript/figures/diagnostics/delayDistBoundary_T' num2str(T) '_N' num2str(N) ...
                '_L' num2str(L) '_a' num2str(alpha,precision) ...
                '_b' num2str(beta,precision)];
            %export_fig([filename '.pdf'],'-depsc')
            set(distributionFig,'PaperUnits','centimeters')
            exportfig(distributionFig,[filename '.eps'],exportOptions);
            system(['epstopdf ' filename '.eps']);
            system(['rm ' filename '.eps']);
        end
        close(distributionFig)
    end
end