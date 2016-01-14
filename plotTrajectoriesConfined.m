% sweep parameters of SPP model and plot individual trajectories
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(2);
betaValues = 2.^(2);
T = 1000;
burnIn = 500;
Trange = burnIn:(burnIn + 100);
NValues = [100 25 9];
LValues = [2 1 0.6];
numAlphas = length(alphaValues);
numBetas = length(betaValues);
bcs = {'noflux', 'noflux', 'free'};

trajectoryFig = figure;
% trajectoryFig.Color='none';
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17,...
    'LineWidth',1);

precision = 2;

%% load results
for permCtr = 1:3
    L = LValues(permCtr);
    N = NValues(permCtr);
    for alphaCtr = 1:numAlphas
        alpha = alphaValues(alphaCtr);
        for betaCtr = 1:numBetas
            beta = betaValues(betaCtr);
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) '_' bcs{1} '-' bcs{2} '-' bcs{3} ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
                '_run1.mat'];
            load(filename)
            % plot trajectories
            plot3(squeeze(cells(:,1,Trange))',squeeze(cells(:,2,Trange))',...
                squeeze(cells(:,3,Trange))','-','LineWidth',1)
            axis equal
            ax = gca;
            ax.XTick = []; ax.YTick = []; ax.ZTick = [];
            ax.XLabel.String = 'x';
            ax.XLabel.Position = ax.XLabel.Position.*[2 0 1];
            ax.YLabel.String = 'y';
            ax.YLabel.Position = ax.YLabel.Position.*[0 2 1];
            ax.ZLabel.String = 'z';
            ax.Box = 'on';
%             ax.Title.String = ['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3)];
%             ax.Title.FontWeight = 'normal';
            % to illustrate confinement, make grey "walls", but leave free
            % boundaries empty (white)
            hold on
            surf(ax.XLim,ax.YLim(2)*[1 1],ax.ZLim'*[1 1],'FaceColor',0.9*[1 1 1])
            surf(ax.XLim(2)*[1 1],ax.YLim,[1 1]'*ax.ZLim,'FaceColor',0.9*[1 1 1])
            % plot a scale bar
            plot3([-0.5 0] + ax.XLim(2),min(ax.YLim)*[1 1],ax.ZLim(1)*[1 1],'k-','LineWidth',3)
            hold off
            %% export figure
            filename = ['manuscript/figures/trajectories_T' num2str(T) '_N' num2str(N) ...
                '_L' num2str(L,'%1.0e') '_' bcs{1} '-' bcs{2} '-' bcs{3} ...
                '_a' num2str(alpha,precision) '_b' num2str(beta,precision)];
            set(trajectoryFig,'PaperUnits','centimeters')
            exportfig(trajectoryFig,[filename '.eps'],exportOptions);
            system(['epstopdf ' filename '.eps']);
            system(['rm ' filename '.eps']);
        end
    end
end
