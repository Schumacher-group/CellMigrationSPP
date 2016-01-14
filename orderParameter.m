function order = orderParameter(cellsTimeSeries)

% input:
% cellsTimeSeries: Array containing the position, and movement direction for
% every cell and time-point. Format is N by [x,y,z] by T.

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

velocityTimeSeries = diff(cellsTimeSeries(:,[x y z],:),1,3);
order = squeeze(sqrt(sum(...
    sum(velocityTimeSeries,1)...% sum velocities of all cells
    .^2,2))...% magnitude of total velocity vector
    ./sum(...
    sqrt(sum(velocityTimeSeries.^2,2))... % magnitude of individual velocities
    ,1)); % sum individual magnitudes over all cells
end