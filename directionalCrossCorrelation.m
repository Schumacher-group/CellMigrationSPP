function crossCorrelation = directionalCrossCorrelation(cellsTimeSeries,lag,r0)

% input:
% cellsTimeSeries: Array containing the position, and movement direction for
% every cell and time-point. Format is N by [x,y,z] by T.

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;

T = size(cellsTimeSeries,3);
N = size(cellsTimeSeries,1);
crossCorrelation = NaN(N*(N-1)/2,T-1-abs(lag));
% only need to go over each pair once, since Cij(tau) = Cji(-tau)

velocityTimeSeries = diff(cellsTimeSeries(:,[x y z],:),1,3);
% normalise for velocity magnitude (should be constant)
directionalTimeSeries = velocityTimeSeries./...
    repmat(sqrt(sum(velocityTimeSeries.^2,2)),1,3);
if lag>= 0
    for timeCtr = 1:(T - 1 - lag)
        crossCorrelation(:,timeCtr) = ...
            squareform(tril(...
            directionalTimeSeries(:,:,timeCtr)*directionalTimeSeries(:,:,timeCtr + lag)'...
            ,-1));
    end
elseif lag<0
    for timeCtr = (1 + abs(lag)):(T - 1)
        crossCorrelation(:,timeCtr - abs(lag)) = ...
            squareform(tril(...
            directionalTimeSeries(:,:,timeCtr)*directionalTimeSeries(:,:,timeCtr + lag)'...
            ,-1));
    end
end