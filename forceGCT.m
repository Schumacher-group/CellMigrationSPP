function force = forceGCT(x,re,ra,r0)
% piecewise force law as in Gregoire et al. Physica D (2003)
force = (x - re)./(ra - re)/4;
force(x>ra) = 1;
force(x>r0) = 0;