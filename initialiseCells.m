function [outArray] = initialiseCells(inArray,L)
% initialises cell positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - initial positions do not currently respect volume exclusion
% - could implement non-random initial positions, e.g. regularly spaced

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3; 
theta = 4; 
phi =   5;

N = size(inArray,1);
if size(L,1)>size(L,2)
    L = L';
end

for cellIdx = 1:N
    % Position, should work for both scalar and vector L
    inArray(cellIdx,[x y z],1) = L.*rand(1,3);
    
    % Direction
    % generate random points on sphere using equal-area projection
    % http://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
    inArray(cellIdx,theta,1) = pi*(2*rand - 1);   % theta between -pi and pi
    inArray(cellIdx,phi,1) = acos(2*rand - 1); % z between -1 and 1
end

outArray = inArray;

