function arrayOut = updateDirection(arrayNow,arrayPrev,L,alpha,beta,...
    r0,ra,re,rc,eta,bc,selfAlign,Nl)
% updates cell directions according to alignment rules and intercellular
% forces

% issues/to-do's:
% - what to use for alignment radius? re? ra? ri?
% - does it matter if we have discontinuities in our force laws curves?
% - could also implement collisions as stopping or random direction
% - periodic boundaries are only implemented for L>2*r0
% - mixed periodic boundary conditions can be quite slow
% - scale Fa by v0 or not?

% short-hand for indexing coordinates
x =     1;
y =     2;
z =     3;
theta = 4;
phi =   5;

N = size(arrayPrev,1);

% find distances between all pairs of cells
ndim = 3;
distanceMatrixXYZ = NaN(N,N,ndim);
for dimCtr = [x y z]
    distanceMatrixXYZ(:,:,dimCtr) = arrayPrev(:,dimCtr*ones(1,N)) - arrayPrev(:,dimCtr*ones(1,N))';
    % reset some distances if boundaries are periodic
    if (~iscell(bc)&&strcmp(bc,'periodic'))||(iscell(bc)&&strcmp(bc{dimCtr},'periodic'))
        if numel(L)==3 % vector domain size [L_x L_y L_z]
            [mirrorIndxRow, mirrorIndxCol] = find(abs(distanceMatrixXYZ(:,:,dimCtr))>=L(dimCtr)/2);
            distanceMatrixXYZ(mirrorIndxRow, mirrorIndxCol,dimCtr) = L(dimCtr) - ...
                distanceMatrixXYZ(mirrorIndxRow, mirrorIndxCol,dimCtr);
        else
            [mirrorIndxRow, mirrorIndxCol] = find(abs(distanceMatrixXYZ(:,:,dimCtr))>=L/2);
            distanceMatrixXYZ(mirrorIndxRow, mirrorIndxCol,dimCtr) = L - ...
                distanceMatrixXYZ(mirrorIndxRow, mirrorIndxCol,dimCtr);
        end
    end
end
distanceMatrix = sqrt(sum(distanceMatrixXYZ.^2,3));

% compute the delaunay triangulation of all cells (to be able to identify
% topological neighbours)
dT3 = delaunayTriangulation(arrayPrev(:,[x y z]));

for cellIdx = 1:N
    % find neighbours in delaunay triangulation
    delaunayNbrs = dT3.isConnected(cellIdx*ones(N,1),(1:N)');
    % calculate force contributions
    
    % alignment force
    % calculate average direction of all cells within re (excluding self unless selfAlignment=true)
    if cellIdx>Nl % standard neighbour alignment for non-informed cells
        alignmentNbrs = delaunayNbrs&distanceMatrix(:,cellIdx)<=re;
        if selfAlign
            alignmentNbrs(cellIdx) = true; % include self-alignment (persistence)
        end
        nbrsPhi = arrayPrev(alignmentNbrs,phi);
        nbrsTheta = arrayPrev(alignmentNbrs,theta);
        Fa = sum([sin(nbrsPhi).*cos(nbrsTheta), sin(nbrsPhi).*sin(nbrsTheta), cos(nbrsPhi)],1)';
    else % the first N cells are informed and align with +x-direction
        Fa = [1; 0; 0];
    end
    % core repulsion (volume exclusion)
    collisionNbrs = distanceMatrix(:,cellIdx)<=rc;
    collisionNbrs(cellIdx) = false; % no self-repulsion
    Nc = nnz(collisionNbrs);
    Fc = sum(... % sum over all collision neighbours
        reshape(distanceMatrixXYZ(cellIdx,collisionNbrs,:),Nc,ndim)... %direction FROM neighbours TO cell (use reshape rather than squeeze for case when Nc = 1)
        ./distanceMatrix(collisionNbrs,cellIdx*ones(1,ndim))... % normalise for distance
        ,1)'; % made into column vector
    
    % intercellular forces
    interactionNbrs = delaunayNbrs&distanceMatrix(:,cellIdx)<=r0&~collisionNbrs;
    interactionNbrs(cellIdx) = false; % no self-interaction
    Ni = nnz(interactionNbrs);
    Fi = (-forceLJS(distanceMatrix(cellIdx,interactionNbrs),re,ra,r0)... % force law, sign-fixed(?)
        *...% sum over all neighbours
        (reshape(distanceMatrixXYZ(cellIdx,interactionNbrs,:),Ni,ndim)... % direction FROM neighbours TO cell (use reshape rather than squeeze for case when Ni = 1)
        ./distanceMatrix(interactionNbrs,cellIdx*ones(1,ndim)))...  % normalise for distence
        )'; % made into column vector
    
    % noise term
    % generate random points on sphere using equal-area projection
    % http://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
    thetaRand = pi*(2*rand - 1);   % theta- between -pi and pi
    zRand = 2*rand - 1; % z between -1 and 1
    % noise scales with number of neighbours (interacting or repulsing),
    % + 1 to include self
    Fnoise = (Ni + Nc + 1)*[sqrt(1 - zRand.^2)*cos(thetaRand);...
        sqrt(1 - zRand.^2)*sin(thetaRand);...
        zRand];
    
    % sum forces
    F = alpha*Fa + beta*(Fi + 1e100*Fc) + eta*Fnoise;
    % 1e100 to represent Inf vector, but still have direction and =0 for beta
    % = 0, should work as long as beta*Fc~beta*O(L)<1.3408e+54
    
    % update directions
    arrayNow(cellIdx,theta) = atan2(F(y),F(x));
    arrayNow(cellIdx,phi) = acos(F(z)/sqrt(sum(F.^2)));
end

arrayOut = arrayNow;