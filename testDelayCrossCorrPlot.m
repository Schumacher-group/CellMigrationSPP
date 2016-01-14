% test delay cross correlation plot
close all
clear

T = 1000;
N = 100;
L = 3;
alpha = 3;
beta = 20;
burnIn = 500;
r0 = 1;
% run simulation
cells = runSPP(T,N,L,alpha,beta,'bc','free');
% discard burn-in
cells = cells(:,:,burnIn:end);

% % plot cell trajectories
plot3(squeeze(cells(:,1,1))',squeeze(cells(:,2,1))',squeeze(cells(:,3,1))',...
    'kx')
hold on
plot3(squeeze(cells(:,1,:))',squeeze(cells(:,2,:))',squeeze(cells(:,3,:))',...
    '-','LineWidth',1)
plot3(squeeze(cells(:,1,end))',squeeze(cells(:,2,end))',squeeze(cells(:,3,end))',...
    'ko')
axis equal
%% analyse trajectories
lagValues = [-100:100];
nLagValues = length(lagValues);
dirCrossCorr = NaN(N*(N-1)/2,nLagValues);
% only need to go over each pair once, since Cij(tau) = Cji(-tau)
minCorr = 0.5;

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
            [pks, locs] = findpeaks(...
                dirCrossCorr(ij,:),lagValues,'NPeaks',1,'SortStr','descend');
            if ~isempty(pks), peakCorrs(ii,jj) = pks; end
            if ~isempty(locs), peakLags(ii,jj) = locs; end
            peakCorrs(jj,ii) = peakCorrs(ii,jj); % correlation is the same for the symmetric pair
            peakLags(jj,ii) = -peakLags(ii,jj);% Cij(tau) = Cji(-tau)
        end
    end
end

figure
% only keep lag times for corr>=minCorr
histogram(peakLags(peakCorrs>=minCorr),20,...
    'DisplayStyle','bar','EdgeColor','none','Normalization','probability');
ax = gca;
ax.Box = 'off';
ax.XLabel.String = '\tau';
ax.YLabel.String = 'P';