% plot an illustrative example of boundary cell detection
close all
clear
%% load data
T = 1000;
burnIn = 500;
tPlot = burnIn + 100;
N = 100;
L = 2;
r0 = 1;
alpha = 4;
beta = 2;
repCtr = 1;
precision = 2;
filename = ['results/' 'T' num2str(T,precision) '_N' num2str(N,precision)...
    '_L' num2str(L,precision) ...
    '_a' num2str(alpha,precision) '_b' num2str(beta,precision) '_run' num2str(repCtr) '.mat'];
out = load(filename);

%% boundary detection as in getBoundaryPoints 
% (here repeated because more output is needed)
maxDist = r0;
points = out.cells(:,:,tPlot);
dt = delaunayTriangulation(points);

% get the connectivity list
clist = dt.ConnectivityList;

% find all edges longer than maxDist
edges = dt.edges;
edgeLengths = sqrt(sum((dt.Points(edges(:,1),:)-dt.Points(edges(:,2),:)).^2,2));
edgesToRemove = edgeLengths>maxDist;
% find triangles/tetrahedra attached to edges
ti = dt.edgeAttachments(edges(edgesToRemove,:));
% remove triangles/tetrahedra
clist(horzcat(ti{:}),:) = [];

% create a new triangulation with the modified connectivity list
nt = triangulation(clist,dt.Points);

%% plot points and new boundary as follows:

plot3(points(:,1),points(:,2),points(:,3),'ko','MarkerSize',3)
hold on
[bt, boundaryPoints] = dt.freeBoundary;
trisurf(bt,boundaryPoints(:,1),boundaryPoints(:,2),boundaryPoints(:,3), ...
    'EdgeColor','r','EdgeAlpha',1,'FaceColor','none')
[bt, boundaryPoints] = nt.freeBoundary;
trisurf(bt,boundaryPoints(:,1),boundaryPoints(:,2),boundaryPoints(:,3), ...
    'EdgeColor',0.3*[1 1 1],'FaceColor','k','FaceAlpha',0.2,'Marker','.','MarkerSize',10)
% format plot
ax = gca;
ax.XTick = []; ax.YTick = []; ax.ZTick = [];
xlabel('x'), ylabel('y'), zlabel('z')
ax.YLabel.Position = ax.YLabel.Position + [0.5, 0.5, 0];
ax.Box = 'on';
ax.ZLim = minmax(points(:,3)');
ax.XLim = minmax(points(:,1)');
ax.YLim = minmax(points(:,2)');
ax.DataAspectRatio = [1 1 1];
% plot a scale bar
plot3([-1 0] + ax.XLim(2),min(ax.YLim)*[1 1],ax.ZLim(1)*[1 1],'k-','LineWidth',3)

%% export figure
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',600,...
    'FontMode','fixed',...
    'FontSize',10);
filename = ['manuscript/figures/boundaryDetection_t' num2str(tPlot) '_N' num2str(N) ...
    '_L' num2str(L,precision) ...
    '_a' num2str(alpha,precision) '_b' num2str(beta,precision)];
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);