% test SPP model
clear
close all


T = 600;
N = 100;
L = 2;
alpha = 4;
beta = 1;
Trange = 500:T;

cells = runSPP(T,N,L,alpha,beta,'bc','free');
order = orderParameter(cells);

plot3(squeeze(cells(:,1,Trange))',squeeze(cells(:,2,Trange))',squeeze(cells(:,3,Trange))','-','LineWidth',1)

axis equal

