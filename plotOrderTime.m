% sweep parameters of SPP model and plot order param over time
close all
clear

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

alphaValues = 2.^(0:7);
betaValues = 2.^(-4:7);
numRepeats = 10;
T = 1000;
burnIn = 500;
N = 100;
L = 2;
numAlphas = length(alphaValues);
numBetas = length(betaValues);

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1);
precision = 2;

%% load results
for alphaCtr = 1:numAlphas
    alpha = alphaValues(alphaCtr);
    for betaCtr = 1:numBetas
        beta = betaValues(betaCtr);
        orderFig = figure;
        orderFig.Color='none';
        hold on
        for repCtr = 1:numRepeats
        % load results
        filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
            '_L' num2str(L,precision) '_a' num2str(alpha,precision) ...
            '_b' num2str(beta,precision) '_run' num2str(repCtr) '.mat'];
        load(filename)
        % plot order parameter over time
        order = orderParameter(cells);
        plot(order)
        end
        ax = gca;
        ax.Title.String = ['\alpha = ' num2str(alpha,3) ', \beta = ' num2str(beta,3)];
        ax.Title.FontWeight = 'normal';
        ax.XLabel.String = 't'; ax.YLabel.String = '\Phi';
        box off
        %% export figure
        filename = ['manuscript/figures/diagnostics/orderTime_T' num2str(T) '_N' num2str(N) ...
            '_L' num2str(L) '_a' num2str(alpha,precision) ...
            '_b' num2str(beta,precision)];
        %         export_fig([filename '.pdf'],'-opengl')
        exportfig(orderFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['rm ' filename '.eps']);
        close(orderFig)
    end
end
