% calculate, plot, and save distribution of peak delay times
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(2:7);
betaValues = 2.^(0:7);
numRepeats = 10;
T = 1000;
burnIn = 500;
N = 100;
Nl = 10;
L = 2;
r0 = 1;
numAlphas = length(alphaValues);
numBetas = length(betaValues);
bcs = 'free';

binWidth = 3;
maxLag = 20;
lagValues = -27:27; % should allow for an integer number of bins
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
histPeakVar = NaN(numAlphas,numBetas,numRepeats);

for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        order = NaN(numRepeats,1);
        for repCtr = 1:numRepeats
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                        '_Nl' num2str(Nl) '_L' num2str(L,precision) '_' bcs ...
                        '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
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
                    histPeakVar(alphaCtr,betaCtr,repCtr) = var(peakLags(threshPeaks));
                end
            end
        end
    end
end
%% plot multi-line diagram of variance
exportOptions.FontSize = 14;
exportOptions.Width = '13';
exportOptions.Height = '11.2';
legendString = num2str(round(alphaValues)');

sigmaFig = figure;
boundedline(betaValues,nanmean(sqrt(histPeakVar),3),...
    permute(nanstd(sqrt(histPeakVar),0,3)./...
    sqrt(sum(~isnan(histPeakVar),3)),[3 2 1]),'alpha','.-','nan','gap')
ax1 = gca;
ax1.XScale = 'log';
ax1.XLabel.String = 'attraction-repulsion strength \beta'; ax1.YLabel.String = 'heterogeneity \sigma(\tau_C)';
absFig.Color='none'; ax1.Box = 'on';
xlim([1 max(betaValues)])
ylim([0 25])
% make an inset of non-log scale
inset = axes('position',...
    [0.28 0.57 0.4 0.345]);
boundedline(betaValues,nanmean(sqrt(histPeakVar),3),...
    permute(nanstd(sqrt(histPeakVar),0,3)./...
    sqrt(sum(~isnan(histPeakVar),3)),[3 2 1]),'alpha','.-','nan','gap')
inset.XLim = [1 max(betaValues)];
inset.YLim = [0 16];
inset.YTick = [0 5 10 15];
inset.Box = 'on';
legH = legend(ax1,legendString,'Location','SouthEast');
legH.Title.String = '\alpha';
legH.Position = [0.7325 0.671 0.1 0.14];
% save figure
filename = ['manuscript/figures/varianceDelayDiagram_T' num2str(T) '_N' num2str(N) ...
    '_Nl' num2str(Nl) '_L' num2str(L)];
set(sigmaFig,'PaperUnits','centimeters')
exportfig(sigmaFig,[filename '.eps'],exportOptions);
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
    '_Nl' num2str(Nl) '_L' num2str(L)];
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);
%%
save(['delayCorrResults' '_Nl' num2str(Nl) '_selfAlign0.mat'],'alphaValues','betaValues','histPeakVar')