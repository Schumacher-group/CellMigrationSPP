% sweep parameters of SPP model and plot individual trajectories next to
% each other
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(4:-1:1);
betaValues = 2.^(-2:2:8);
T = 1000;
burnIn = 500;
Trange = burnIn:(burnIn + 100);
N = 10;
L = [1 1 0.05];
numAlphas = length(alphaValues);
numBetas = length(betaValues);

trajectoryFig = figure;
trajectoryFig.Color='none';
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',20,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1);

bcs = {'free', 'free', 'noflux'};

precision = 2;

%% load results
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        % load results
        filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
            '_L' num2str(L(1),precision) '_' num2str(L(2),precision)...
            '_' num2str(L(3),precision) '_' bcs{1} '-' bcs{2} '-' bcs{3}...
            '_a' num2str(alpha,precision) '_b' num2str(beta,precision) '_run1.mat'];
        load(filename)
        % plot trajectories
        subplot(numAlphas,numBetas,(alphaCtr - 1)*numBetas + betaCtr)
        plot(squeeze(cells(:,1,Trange))',squeeze(cells(:,2,Trange))',...
            '-','LineWidth',0.5)
        axis image
        box on
        ax = gca;
        ax.XTick = []; ax.YTick = []; ax.ZTick = [];
        ax.Title.String = ['\alpha=' num2str(alpha,3) ', \beta=' num2str(beta,3)];
        ax.Title.FontWeight = 'normal';
        %         ax.Visible = 'off';
        if alpha==2&&beta==0.25
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
        end
    end
end
%% export figure
filename = ['manuscript/figures/trajectories_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L(1)) '_semi2D'];
%         export_fig([filename '.pdf'])
exportfig(trajectoryFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);