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
N = 100;
L = 2;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

trajectoryFig = figure;
trajectoryFig.Color='none';
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',17,...
    'LineWidth',1);

precision = 2;

%% load results
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        % load results
        filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
            '_L' num2str(L,precision) '_a' num2str(alpha,precision) ...
            '_b' num2str(beta,precision) ...%'_selfAlign'...
            '_run1.mat'];
        load(filename)
        % plot trajectories
        plot3(squeeze(cells(:,1,Trange))',squeeze(cells(:,2,Trange))',...
            squeeze(cells(:,3,Trange))','-','LineWidth',1)
        axis equal
        ax = gca;
        ax.XTick = []; ax.YTick = []; ax.ZTick = [];
        ax.XLabel.String = 'x';
        ax.XLabel.Position = ax.XLabel.Position.*[1.2 0.95 1];
        ax.YLabel.String = 'y';
        ax.YLabel.Position = ax.YLabel.Position.*[1.1 0.9 1];
        ax.ZLabel.String = 'z';
        %         ax.Title.String = ['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3)];
        %         ax.Title.FontWeight = 'normal';
        % plot a scale bar
        hold on
        plot3([-0.5 0] + ax.XLim(2),min(ax.YLim)*[1 1],ax.ZLim(1)*[1 1],'k-','LineWidth',3)
        hold off
        ax.Box = 'on';
        %% export figure
        filename = ['manuscript/figures/trajectories_T' num2str(T) '_N' num2str(N) ...
            '_L' num2str(L) '_a' num2str(alpha,precision) ...%'_selfAlign' ...
            '_b' num2str(beta,precision)];
        set(trajectoryFig,'PaperUnits','centimeters')
        exportfig(trajectoryFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['rm ' filename '.eps']);
    end
end
