function [ reflected ] = reflectDirection(angles,dimension)
% reflects movement direction upon incidence onto noflux boundary
% angles have to be given in [theta, phi], where theta and phi are column
% vectors
x =     1;
y =     2;
z =     3;

switch dimension
    case x % theta goes to pi - theta (or -pi + theta)
        angles(:,1) = wrapToPi(pi - angles(:,1));
    case y % theta goes to -theta
        angles(:,1) = - angles(:,1);
    case z % phi goes to pi - phi
        angles(:,2) = pi - angles(:,2);
end
reflected = angles;
end

