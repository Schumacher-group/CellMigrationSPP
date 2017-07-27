% A function to run SPP model simulations
% LJ Schumacher January 2016
% With thanks to Yasha Sharma and Diego Vargas for insight into an earlier
% implementation
%
% inputs
% T: simulation duration (number of time-steps)
% N: number of cells
% L: size of region containing cells' initial positions, in units of
%   interaction radii. Can be single number for cubic region, or
%   3 element vector [L_x L_y L_z], should be >2*r0
% alpha: relative weight of alignment interaction
% beta: relative weight of intercellular forces
%
% optional inputs, with default values, partly as in Gregoire et al. Physica D (2003)
% v0: speed (default 0.05)
% r0: interaction radius (default 1.0), should be >ra,re,rc
% ra: attraction radius (default 0.8)
% re: preferred radius (default 0.5)
% rc: core repulsion radius (default 0.2), should be <re,ra,r0
% eta: noise strength (default 1.0)
% bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 3 element array {'bcx','bcy','bcz'} for different
%   bcs along different dimensions
% selfAlign: whether cells should include their own direction when
%  averaging the direction of neighbours. Can be seen as an implementation
%  of persistence. If selfAlign=true, the balance of alpha and eta control 
%  the degree of persistence (default 'false')
% Nl: number of informed/leader cells, which instead of aligning with the
% group, always align with the +x-direction
%
% outputs
% cells: Array containing the position, and movement direction for
% every cell and time-point. Format is N by [x,y,z,theta,phi] by T.

% issues/to-do's:
% - implement checks on relationships btw parameters, e.g. L>2r0, r0>r*,
% rc<r*, ... for noflux min L should be v0*dT to prevent overshoot 

function cells = runSPP(T,N,L,alpha,beta,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) x>0&&~mod(x,1);
checkInt0 = @(x) x>=0&&~mod(x,1);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'L',@checkL);
addRequired(iP,'alpha',@isnumeric);
addRequired(iP,'beta',@isnumeric);
addOptional(iP,'v0',0.05,@isnumeric)
addOptional(iP,'r0',1.0,@isnumeric)
addOptional(iP,'ra',0.8,@isnumeric)
addOptional(iP,'re',0.5,@isnumeric)
addOptional(iP,'rc',0.2,@isnumeric)
addOptional(iP,'eta',1.0,@isnumeric)
addOptional(iP,'bc','free',@checkBcs)
addOptional(iP,'selfAlign','false',@islogical)
addOptional(iP,'Nl',0,checkInt0)
parse(iP,T,N,L,alpha,beta,varargin{:})
v0 = iP.Results.v0;
r0 = iP.Results.r0;
ra = iP.Results.ra;
re = iP.Results.re;
rc = iP.Results.rc;
eta = iP.Results.eta;
bc = iP.Results.bc;
selfAlign = iP.Results.selfAlign;
Nl = iP.Results.Nl;

cells = NaN(N,5,T);
cells = initialiseCells(cells,L);

for t=2:T
    % update direction
    cells(:,:,t) = updateDirection(cells(:,:,t),cells(:,:,t-1),L,alpha,beta,...
    r0,ra,re,rc,eta,bc,selfAlign,Nl);
    % update position
    cells(:,:,t) = updatePosition(cells(:,:,t),cells(:,:,t-1),v0,bc,L);
end
end

function LCheck = checkL(x)
LCheck = false;
if numel(x)==3||numel(x)==1
    LCheck = isnumeric(x);
end
end

function BcCheck = checkBcs(x)
BcCheck = false;
validBcs = {'free','noflux','periodic'};
if iscell(x)&&numel(x)==3
    BcCheck = any(validatestring(x{1},validBcs))...
        &any(validatestring(x{2},validBcs))...
        &any(validatestring(x{3},validBcs));
else
    BcCheck = any(validatestring(x,validBcs));
end
end