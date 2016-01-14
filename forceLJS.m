function force = forceLJS(x,re,ra,r0)
% modified force law to more realistically represent CAM unbinding
force = (x - re)./(ra - re);
force(x>ra) = exp(-2*(x(x>ra)-ra)./(r0 - ra));
force(x>r0) = 0;