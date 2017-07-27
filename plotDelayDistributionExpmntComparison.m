% plot distribution of peak delay times from SPP model for comparison with
% published data
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 4:8;
betaValues = 2.^(1:0.5:6);
numRepeats = 10;
T = 1000;
burnIn = 500;
NValues = [10 20];
LValues = [0.9 1.1];
nTimePoints = 50;
r0 = 1;
numAlphas = length(alphaValues);
numBetas = length(betaValues);
bcs = 'free';

binWidth = 1;
smoothWidth = 3;
maxLag = 20;
lagValues = -22:22; % should allow for an integer number of bins
% use two more bins than wanted (one for nicer plotting, one to discard edges)
nLagValues = length(lagValues);
minCorr = 0.5;
minOrder = 0.79; % choose range representative of experiments
maxOrder = 0.97;

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',12,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',16,...
    'LineWidth',2,...
    'Renderer','opengl');

precision = 2;

experimentalReference = {'delayDistSharmeEtAlG1V2C1_0320_0930.mat',...
    'delayDistSharmeEtAlG2V1C1_0130_0550.mat'};

%% load results
binCentres = (min(lagValues)+2*binWidth):binWidth:(max(lagValues)-2*binWidth);
for permCtr = 1:2
    N = NValues(permCtr);
    L = LValues(permCtr);
    load(experimentalReference{permCtr},'expmntlLagDist');
    % normalise reference distribution
    expmntlLagDist(2,:) = expmntlLagDist(2,:)/sum(expmntlLagDist(2,:));
    bestFit = Inf; bestAlpha = NaN; bestBeta = NaN;
    histPeakLag = NaN(numAlphas,numBetas,numRepeats,(nLagValues - 1)/binWidth - 1);
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
                    '_L' num2str(L,precision) '_' bcs ...
                    '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
                    '_selfAlign0_run' num2str(repCtr) '.mat'];
                out = load(filename);
                % discard burn-in
                out.cells = out.cells(:,:,burnIn:end);
                % choose a time segment of similar length as
                % the experimental time series
                startSample=1;
                out.cells = out.cells(:,:,startSample-1+(1:nTimePoints));
                order(repCtr) = mean(orderParameter(out.cells));
                if order(repCtr) >= minOrder && order(repCtr) <= maxOrder
                    dirCrossCorr = NaN(N*(N-1)/2,nLagValues);
                    % only need to go over each pair once, since Cij(tau) = Cji(-tau)
                    
                    % calculate cross-correlations
                    for lagCtr = 1:nLagValues % can be parfor-ed
                        lag = lagValues(lagCtr);
                        dirCrossCorr(:,lagCtr) = mean(directionalCrossCorrelation(out.cells,lag,r0),2);
                    end
                    
                    % only keep cross-correlations if cells have been neighbours (within r0)
                    pairDistances = NaN(N*(N-1)/2,nTimePoints);
                    for timeCtr=1:nTimePoints % can be parfor-ed
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
                        histPeakLag(alphaCtr,betaCtr,repCtr,:) = ...
                            smooth(histcounts(peakLags(threshPeaks),...
                            'BinWidth',binWidth,'BinLimits',[min(lagValues)+binWidth/2 max(lagValues)-binWidth/2],...
                            'Normalization','probability'),smoothWidth);
                    end
                end
            end
            if plotCurrent
                boundedline(binCentres,...
                    squeeze(nanmean(histPeakLag(alphaCtr,betaCtr,:,2:end-1),3)),...
                    2*squeeze(nanstd(histPeakLag(alphaCtr,betaCtr,:,2:end-1),0,3)),'.-');
                % plot comparison with experimental distributions
                plot(expmntlLagDist(1,:)/10,expmntlLagDist(2,:),'k-')
                % calculate qualitiy of match
                currentFit = mean(sum((squeeze(histPeakLag(alphaCtr,betaCtr,:,2:end-1))...
                    - expmntlLagDist(2*ones(numRepeats,1),:)).^2,2));
                if currentFit<bestFit
                    bestFit = currentFit;
                    bestAlpha = alpha;
                    bestBeta = beta;
                end
                % format figure
                ylim([0 0.16])
                xlim([-maxLag maxLag])
                ax = gca;
                ax.YTick = 0:0.04:0.16;
                ax.YTickLabel = num2str(ax.YTick','%1.2f');
                ax.Box = 'on';
                ax.XLabel.String = 'peak delay time \tau_C';
                ax.YLabel.String = 'probability P';
                ax.Title.String = {['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3),...
                    ', \langle\Phi\rangle = ' num2str(nanmean(order),2) '±' num2str(2*nanstd(order),1)];''};
                ax.Title.FontWeight = 'normal';
                %% export figure
                filename = ['manuscript/figures/diagnostics/delayDist_T' num2str(T) '_N' num2str(N) ...
                    '_L' num2str(L,'%2.0e\n') '_a' num2str(alpha,precision) ...
                    '_b' num2str(beta,precision) '_selfAlign0'];
                %export_fig([filename '.pdf'],'-depsc')
                set(distributionFig,'PaperUnits','centimeters')
                exportfig(distributionFig,[filename '.eps'],exportOptions);
                system(['epstopdf ' filename '.eps']);
                system(['rm ' filename '.eps']);
            end
            close(distributionFig)
        end
    end
    disp(['Best fit for ' experimentalReference{permCtr} ' of ' num2str(bestFit) ...
        ' was achieved with alpha=' num2str(bestAlpha) ', beta=' num2str(bestBeta)])
    %% save results
    save(['delayCorrResults_N' num2str(N) '_L' num2str(L) '_selfAlign0' ...
        '.mat'],'alphaValues','betaValues','histPeakLag','bestAlpha','bestBeta')
end
