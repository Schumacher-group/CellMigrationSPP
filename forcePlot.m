% plot intercellular forces
clear
close all

rc = 0.2;
re = 0.5;
ra = 0.8;
r0 = 1;
dr = 0.01;

f = @(x) forceGCT(x,re,ra,r0);
f2 = 1;
f3 = @(x) forceLJS(x,re,ra,r0);

r1 = rc:dr:ra;
plot(r1,f(r1),'k--')
hold on
plot(r1,f3(r1),'k')

r2 = (ra+eps):dr:r0;
plot(r2,f(r2),'k--')
plot(r2,f3(r2),'k')

plot([rc rc],[-1,f(rc)],'k--')

set(gca,'ytick',[-1,f(rc),0,f(ra),f(r0)],'XAxisLocation','origin')
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.2f'));
set(gca,'xtick',[0,rc,re,ra,r0],'xticklabel',{num2str(0);'r_c';'r_e';'r_a';'r_0'})

xlim([0 r0])
ylabel('force (arbitraty units)')
xlabel('distance (between cells)')
% legend('GCT','SMB','Location','NorthWest')
box off

%% export figure
set(gcf,'Color','none')
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = 'manuscript/figures/forcePlot';
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);