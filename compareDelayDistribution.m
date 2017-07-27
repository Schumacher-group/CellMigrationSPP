 % compare the distributions of peak delay times for different parameters
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^4;
betaValues = 2.^[5 6 7];
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

binCentres = (min(lagValues)+2*binWidth):binWidth:(max(lagValues)-2*binWidth);
plotColors = lines(numBetas*numAlphas);
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',9,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

precision = 2;

%% load results
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    histPeakLag = NaN(numBetas,numRepeats,(nLagValues - 1)/binWidth - 1);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        for repCtr = 1:numRepeats
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) '_a' num2str(alpha,precision) ...
                '_b' num2str(beta,precision) '_run' num2str(repCtr) '.mat'];
            load(filename)
            % discard burn-in
            cells = cells(:,:,burnIn:end);
            order = mean(orderParameter(cells));
            if order >= minOrder
                dirCrossCorr = NaN(N*(N-1)/2,nLagValues);
                % only need to go over each pair once, since Cij(tau) = Cji(-tau)
                
                % calculate cross-correlations
                for lagCtr = 1:nLagValues
                    lag = lagValues(lagCtr);
                    dirCrossCorr(:,lagCtr) = mean(directionalCrossCorrelation(cells,lag,r0),2);
                end
                
                % only keep cross-correlations if cells have been neighbours (within r0)
                pairDistances = NaN(N*(N-1)/2,T-burnIn);
                for timeCtr=1:T-burnIn
                    pairDistances(:,timeCtr) = pdist(cells(:,1:3,timeCtr));
                end
                notNeighbours = min(pairDistances,[],2)>r0;
                dirCrossCorr(notNeighbours,:) = NaN;
                
                % find and plot peaks in crosscorr curves
                peakLags = NaN(N);
                peakCorrs = NaN(N);
                
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
                % warnigns                
                threshPeaks = peakCorrs>=minCorr;
                % discard peak lag times outside relevant range from
                % further calculation
                threshPeaks(abs(peakLags(threshPeaks))>maxLag) = false;
                if any(threshPeaks(:))
                    histPeakLag(betaCtr,repCtr,:) = histcounts(peakLags(threshPeaks),...
                        'BinWidth',binWidth,'BinLimits',[min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                        'Normalization','probability');
                end
            end
        end
    end
    %% plot mean distribution and bounds
    distributionsFig = figure;
    hold on
    boundedline(binCentres,squeeze(nanmean(histPeakLag(:,:,2:end-1),2)),...
        permute(nanstd(histPeakLag(:,:,2:end-1),0,2),[3 2 1]),'alpha','.-');
    ax = gca;
    ax.Box = 'on';
    ax.XLabel.String = 'peak delay time \tau_C';
    ax.YLabel.String = 'probability P';
    % adjust ticks and significant figures
    ax.YTick = 0:0.1:03;
    ax.YTickLabel = num2str(ax.YTick','%1.1f');
    title(['\alpha = ' num2str(alpha,3)],'FontWeight','normal')
    legH = legend(num2str(betaValues',3));
    legH.Title.String = '\beta';
    ylim([0 max(max(max(histPeakLag)))])
    xlim([-20 20])
    %% export figure
    set(distributionsFig,'PaperUnits','centimeters')
    filename = ['manuscript/figures/delayDist_T' num2str(T) '_N' num2str(N) ...
        '_L' num2str(L) '_a' num2str(alpha,precision)];
    exportfig(distributionsFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']);
    close(distributionsFig)
end
