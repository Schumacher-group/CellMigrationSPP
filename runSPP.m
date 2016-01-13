% A function to run SPP model simulations
% LJ Schumacher January 2016
% With thanks to Yasha Sharma and Diego Vargas for insight into an earlier
% implementation

% inputs
% T: simulation duration (number of time-steps)
% N: number of cells
% L: size of region containing cells' initial positions, in units of
%   interaction radii. Can be single number for cubic region, or
%   3 element vector [L_x L_y L_z]
% alpha: relative weight of alignment interaction
% beta: relative weight of intercellular forces
%
% optional inputs, mostly as in Gregoire et al. Physica D (2003)
% v0: speed (default 0.05)
% r0: interaction radius (default 1.0)
% ra: attraction radius (default 0.8)
% re: preferred radius (default 0.5)
% rc: core repulsion radius (default 0.2)
% eta: noise strength (default 1.0)
% bc: boundary conidition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 3 element array {'bcx','bcy','bcz'} for different
%   bcs along different dimensions
%
% outputs
% cells: Matrix containing the index, position, and movement direction for
% every cell and time-point. Format is ...

function cells = runSPP(T,N,L,alpha,beta,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) x>0&&~mod(x,1);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'L',@isnumeric);
addRequired(iP,'alpha',@isnumeric);
addRequired(iP,'beta',@isnumeric);
addOptional(iP,'v0',0.05,@isnumeric)
addOptional(iP,'r0',1.0,@isnumeric)
addOptional(iP,'ra',0.8,@isnumeric)
addOptional(iP,'re',0.5,@isnumeric)
addOptional(iP,'rc',0.2,@isnumeric)
addOptional(iP,'eta',1.0,@isnumeric)
addOptional(iP,'bc','free',@checkBcs)
parse(iP,T,N,L,alpha,beta,varargin{:})

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
