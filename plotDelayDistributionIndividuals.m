% sweep parameters of SPP model and plot distribution of peak delay times
% for individual cells
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 16%2.^(2:1:7);
betaValues = 32%2.^(1:1:7);
numRepeats = 10;
T = 1000;
burnIn = 500;
N = 100;
L = 2;
r0 = 1;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

binWidth = 1;
smoothWidth = 3;
maxLag = 20;
lagValues = -24:24; % should allow for an integer number of bins
% use two more bins than wanted (one fpr nicer plotting, one to discard edges)
nLagValues = length(lagValues);
minCorr = 0.5;
minOrder = 0.1; % as computing correlations for all pairs and delays is expensive, only do for relatively ordered collectives

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',30.5,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1);

precision = 2;

%% load results
binCentres = (min(lagValues)+2*binWidth):binWidth:(max(lagValues)-2*binWidth);
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        distributionFig = figure;
        distributionFig.Color='none';
        for repCtr = 1:(numRepeats-1)
            subplot(3,3,repCtr), hold on
            plotCurrent = 0;
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ... %'_selfAlign' ...
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
                % %                 boundaryMtx = false(N,T-burnIn);
                for timeCtr=1:T-burnIn % can be parfor-ed
                    pairDistances(:,timeCtr) = pdist(out.cells(:,1:3,timeCtr));
                    % %                     % check which cells are at the boundary
                    % %                     [~, boundaryPointIDs] = getBoundaryPoints(out.cells(:,1:3,timeCtr),r0);
                    % %                     boundaryMtx(boundaryPointIDs,timeCtr) = true;
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
                    [thrPkrow, thrPkcol] = find(threshPeaks);
                    % calculate all individual delay distributions
                    delayDist = NaN(N,length(binCentres) + 2);
                    asymmetries = NaN(N,1);
                    for cellIdx=1:N
                        cellIdcs = thrPkcol==cellIdx;
                        cellLags = peakLags(thrPkrow(cellIdcs),thrPkcol(cellIdcs));
                        if nnz(~isnan(cellLags))>1
                            delayDist(cellIdx,:) = smooth(histcounts(cellLags(:),...
                                'BinWidth',binWidth,'BinLimits',...
                                [min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                                'Normalization','probability'),smoothWidth);
                            % calculate asymmetry (away from zero, not it's mean)
                            asymmetries(cellIdx) = min(1,abs(binCentres*delayDist(cellIdx,2:end-1)')/maxLag);
                        end
                    end
                    if any(any(delayDist))
                        % plot distributions in increasing order of asymmetry
                        [~, plotOrder] = sort(asymmetries);
                        for cellIdx=plotOrder'
                            if ~all(isnan(delayDist(cellIdx,:)))
                                % set color and transparency
                                plotColor = [0.33 0.33 0.33 sqrt(asymmetries(cellIdx))]*(1 - nthroot(asymmetries(cellIdx),3)) ...
                                    + [1 0 0 sqrt(asymmetries(cellIdx))]*nthroot(asymmetries(cellIdx),3);
                                plot(binCentres,delayDist(cellIdx,2:end-1),'Color',plotColor)
                            end
                        end
                        % format figure
                        ax = gca;
                        ax.Box = 'on';
                        ax.XLim = [-maxLag maxLag];
                        ax.XTick = [-maxLag:5:maxLag];
%                         ax.XTickLabel{2} = '\tau_C';
                        ax.XLabel.String = 'peak delay time \tau_C';
                        ax.YLim = [0 ceil(max(max(delayDist(:,2:end-1)))*10)/10];
                        ax.YTick = linspace(0,ax.YLim(2),3);
                        ax.YTickLabel = num2str(ax.YTick','%1.1f');
%                         ax.YTickLabel{2} = 'P';
                        ax.YLabel.String = 'probability P';
                        ax.Title.String = ['individual distributions, \alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3)];
                        ax.Title.FontWeight = 'normal';
                    end
                end
            end
        end
        
        %% export figure
        filename = ['manuscript/figures/diagnostics/delayDistIndividuals_T' num2str(T) '_N' num2str(N) ...
            '_L' num2str(L) '_a' num2str(alpha,precision) ...
            '_b' num2str(beta,precision) ... %'_selfAlign'
            ];
        %export_fig([filename '.pdf'],'-depsc')
        set(distributionFig,'PaperUnits','centimeters')
        exportfig(distributionFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['rm ' filename '.eps']);
        close(distributionFig)
    end
end