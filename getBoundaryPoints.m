function [boundaryPoints, boundaryPointIDs] = getBoundaryPoints(points,maxDist)
% gives the points on the boundary of a cloud of points, as determined
% by a delaunay triangulation, from which we remove any edges that are
% longer than maxDist (which will not be a delauny triangulation
% anymore!)

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

% determinte the free boundary of the new triangulation
[~, boundaryPoints] = nt.freeBoundary;
boundaryPointIDs = find(all(ismember(points,boundaryPoints),2));
% plot points and new boundary as follows:
% plot3(points(:,1),points(:,2),points(:,3),'ko','MarkerSize',3)
% hold on
% [bt, boundaryPoints] = dt.freeBoundary;
% trisurf(bt,boundaryPoints(:,1),boundaryPoints(:,2),boundaryPoints(:,3), ...
%     'EdgeColor','r','EdgeAlpha',1,'FaceColor','none')
% [bt, boundaryPoints] = nt.freeBoundary;
% trisurf(bt,boundaryPoints(:,1),boundaryPoints(:,2),boundaryPoints(:,3), ...
%     'EdgeColor','0.3*[1 1 1],'FaceColor','k','FaceAlpha',0.2,'Marker','.','MarkerSize',10)
end