% calculate and plot order phase diagrams
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(0:0.5:7);
betaValues = 2.^(-4:0.5:7);
numRepeats = 10;
T = 1000;
burnIn = 500;
N = 100;
L = 2;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

precision = 2;
order = NaN(numAlphas, numBetas, numRepeats);
orderStd = NaN(size(order));
%% load results
for repCtr = 1:numRepeats
    for alphaCtr = 1:numAlphas
        alpha = alphaValues(alphaCtr);
        for betaCtr = 1:numBetas
            beta = betaValues(betaCtr);
            % load results
            filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
                '_L' num2str(L,precision) '_a' num2str(alpha,precision) ...
                '_b' num2str(beta,precision) '_run' num2str(repCtr) '.mat'];
            load(filename)
            % calculate order
            order(alphaCtr,betaCtr,repCtr) = mean(orderParameter(cells(:,:,burnIn:end)));
            orderStd(alphaCtr,betaCtr,repCtr) = std(orderParameter(cells(:,:,burnIn:end)));
        end
    end
end

%% plot phase diagram 
orderFig = figure;
nContours = 10;
contourf(betaValues,alphaValues,mean(order,3),nContours,'LineColor','none')
ax1 = gca;
ax1.XScale = 'log'; ax1.YScale = 'log';
ax1.XLabel.String = 'attraction-repulsion strength \beta'; ax1.YLabel.String = 'alignment strength \alpha';
ax1.XTick = [0.1 1 10 100];
ax1.XLim = 2.^[-1 7];
ax1.DataAspectRatio = [1 1 1];
cbar1 = colorbar;
cbar1.Label.String = 'directional order \Phi';
cbar1.YTick = [0.0:0.2:1.0];
% reformat significant figures
cbar1.YTickLabel = num2str(cbar1.YTick','%1.1f');
colormap(flipud(gray(nContours)))
caxis([0 1])

meanErrFigure = figure;
contourf(betaValues,alphaValues,std(order,[],3)./sqrt(numRepeats),nContours,'LineColor','none')
ax2 = gca;
ax2.XScale = 'log'; ax2.YScale = 'log';
ax2.XLabel.String = '\beta'; ax2.YLabel.String = '\alpha';
cbar2 = colorbar;
cbar2.Label.String = '\Phi';
colormap(cool(nContours))

%% export figure
orderFig.Color='none';
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',12,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',9,...
    'LineWidth',2,...
    'Renderer','opengl');
filename = ['manuscript/figures/orderDiagram_T' num2str(T) '_N' num2str(N) ...
    '_L' num2str(L)];
set(orderFig,'PaperUnits','centimeters')
exportfig(orderFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);