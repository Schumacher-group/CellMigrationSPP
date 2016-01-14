% sweep parameters of SPP model and plot individual trajectories next to
% each other
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(3:-1:1);
betaValues = 2.^(0:5);
T = 1000;
burnIn = 500;
Trange = burnIn:(burnIn + 100);
N = 100;
Nl = 10;
L = 2;
numAlphas = length(alphaValues);
numBetas = length(betaValues);
bcs = 'free';
selfAlign = false;

trajectoryFig = figure;
trajectoryFig.Color='none';
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',20,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',9,...
    'LineWidth',1);

precision = 2;

%% load results
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        % load results
        filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                        '_Nl' num2str(Nl) '_L' num2str(L,precision) '_' bcs ...
                        '_a' num2str(alpha,precision) '_b' num2str(beta,precision) ...
                        '_selfAlign' num2str(selfAlign) '_run1.mat'];
        load(filename)
        % plot trajectories
        subplot(numAlphas,numBetas,(alphaCtr - 1)*numBetas + betaCtr)
        plot3(squeeze(cells(:,1,Trange))',squeeze(cells(:,2,Trange))',...
            squeeze(cells(:,3,Trange))','-','LineWidth',0.5)
        axis image
        box on
        ax = gca;
        ax.XTick = []; ax.YTick = []; ax.ZTick = [];
        ax.Title.String = {['\alpha=' num2str(alpha,3) ', \beta=' num2str(beta,3)];...
            %['\Phi=' num2str(mean(orderParameter(cells(:,:,Trange))),1)]...
            };
        ax.Title.FontWeight = 'normal';
        ax.Position= ax.Position.*[1 1 1.2 1.2]; % stretch its width and height to reduce whitespace btw subplots
        if alpha==1&&beta==1
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
        end
        %         ax.Visible = 'off';
        % plot a scale bar
        hold on
        plot3([-1 0] + ax.XLim(2),min(ax.YLim)*[1 1],ax.ZLim(1)*[1 1],'k-','LineWidth',2)
        hold off
    end
end
%% export figure
filename = ['manuscript/figures/trajectories_T' num2str(T) '_N' num2str(N) ...
    '_Nl' num2str(Nl) '_L' num2str(L) '_' bcs '_selfAlign' num2str(selfAlign)];
set(trajectoryFig,'PaperUnits','centimeters')
exportfig(trajectoryFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);

system(['rm ' filename '.eps']);