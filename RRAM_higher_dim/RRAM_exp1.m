% This script generates the first experiment in the manuscript with the
% RRAM (section 4.3.1).

clear

d = 4;
nn = 20;
n = nn*ones(1,d);

rA = [1 4*ones(1,d-1) 1];
rng(1)
A = TTeMPS_randn(rA, n);
A_full = full(A);

rx = [1, ones(1,d-1), 1];
X0 = TTeMPS_randn(rx, n);

sizeOmega = round(0.1*prod(n));
sizeGamma = round(sizeOmega/4);
rng(1)
Omega_Gamma_ind = randperm(prod(n),sizeOmega+sizeGamma)';
Omega_Gamma = ind2sub2(n,Omega_Gamma_ind);
Omega = Omega_Gamma(1:sizeOmega,:);
Gamma = Omega_Gamma(sizeOmega+1:end,:);
Omega_ind = sub2ind2(n,Omega);
Gamma_ind = sub2ind2(n,Gamma);
sizeOmega_C = []; Omega_C = {[],[],[],[]}; Omega_C_ind = []; A_Omega_C = [];
if isnumeric(A)
    A_Omega = A(Omega_ind);
    A_Gamma = A(Gamma_ind);
else
    A_Omega = A(Omega);
    A_Gamma = A(Gamma);
end

%%
opts1 = struct('maxiter',15,'maxiter_final',15,'gradtol',10^(-8),...
    'reltol',10^(-8),'tol',10^-8);
opts1.maxrank = max(rA);
opts2 = struct('kmax',15,'eps_gamma',1,...
    'r_max',10*ones(1,d-1),'s_max',8*ones(1,d-1),'Delta',0.8);

opts2.maxiter = opts1.maxiter;
opts2.gradtol = opts1.gradtol;
opts2.eps_omega = opts1.tol;
opts2.tol = opts1.tol;
opts2.reltol = opts1.reltol;

%%
[X2,cost_tmp2,test_tmp2,stats2,ranks2] = RRAM_TT_completion_gen(...
    A_Omega, Omega, A_Gamma, Gamma, X0, opts2 );

%%
[X,cost_tmp,test_tmp,stats1,ranks] = completion_rankincrease(...
    'GeomCG', A_Omega, Omega, A_Gamma, Gamma, A_Gamma, Gamma, X0, opts1); 

%%
figure
t = tiledlayout(3,5, 'Padding', 'none', 'TileSpacing', 'compact'); 
stats1_r = cumsum(stats1.rankidx); stats2_r = cumsum(stats2.rankidx);

set(0,'defaultAxesFontSize',24)
nexttile([3 2])
semilogy(1:length(cost_tmp),0.5*cost_tmp.^2,'b-+','Markersize',8,'Linewidth',2);
hold on
semilogy(1:length(cost_tmp2),0.5*cost_tmp2.^2,'r-o','Markersize',8,'Linewidth',2);
axis tight
title('$f_{\Omega}\big(X^{(j)}\big)/||A_{\Omega}||^2$','interpreter','latex')

y = ylim(gca);
if y(1) < 10^-1
    yticks(10.^(ceil(log10(y(1))):2:0))
end
for j=1:length(stats1_r)
    semilogy([stats1_r(j), stats1_r(j)],y,'--','color',[0.7,0.7,0.7]);
end
for j=1:length(stats2_r)
    semilogy([stats2_r(j), stats2_r(j)],y,'--','color',[0.7,0.7,0.7]);
end

xlabel('$j$, inner (cumulative)','interpreter','latex') 

%%   
nexttile([3 2])
p1 = semilogy(stats1.time,0.5*test_tmp.^2,'b-+','Markersize',8,'Linewidth',2);
hold on
p2 = semilogy(stats2.time,0.5*test_tmp2.^2,'r-o','Markersize',8,'Linewidth',2);

L{1} = '$\mathrm{RRAM}_{\mathrm{randn}}$';
L{2} = '$\mathrm{RRAM}_{\tilde{\mathcal{P}}}$';
axis tight

y = ylim(gca); x = xlim(gca);
if y(1) < 0.1
    yticks(10.^(ceil(log10(y(1))):2:0))
end
for j=1:length(stats1.time)
    semilogy( [stats1.time(j), stats1.time(j)], y, '--','color',[0.7,0.7,0.7]);
end
for j=1:length(stats2.time)
    semilogy( [stats2.time(j), stats2.time(j)], y, '--','color',[0.7,0.7,0.7]);
end

xticks([0:2:x(2)])
axis tight 
legend([p1,p2],L,'interpreter','latex','location','southeast')    
xlabel('seconds','interpreter','latex')
title('$f_{\Gamma}\big(X^{(k)}\big)/||A_{\Gamma}||^2$','interpreter','latex')


%%
p = nexttile;
plot(ranks(:,2),'b-+','Markersize',8,'Linewidth',2);
hold on
plot(ranks2(:,2),'r-o','Markersize',8,'Linewidth',2);
axis tight

y = ylim(p); x = xlim(p);
for j=2:x(2)-1
    semilogy([j,j], y, '--','color',[0.7,0.7,0.7]);
end

xlabel('$k$, outer','interpreter','latex')
ylabel('$r_1$','interpreter','latex')
yticks([0:2:y(2)])

%%    
p = nexttile;
plot(ranks(:,3),'b-+','Markersize',8,'Linewidth',2);
hold on
plot(ranks2(:,3),'r-o','Markersize',8,'Linewidth',2);
axis tight

y = ylim(p);
x = xlim(p);
for j=2:x(2)-1
    semilogy([j,j], y, '--','color',[0.7,0.7,0.7]);
end

yticks([0:2:y(2)])
xlabel('$k$, outer','interpreter','latex')
ylabel('$r_2$','interpreter','latex')

%%
p = nexttile;
plot(ranks(:,4),'b-+','Markersize',8,'Linewidth',2);
hold on
plot(ranks2(:,4),'r-o','Markersize',8,'Linewidth',2);
axis tight

y = ylim(p);
x = xlim(p);
for j=2:x(2)-1
    semilogy([j,j], y, '--','color',[0.7,0.7,0.7]);
end

yticks([0:2:y(2)])
xlabel('$k$, outer','interpreter','latex')
ylabel('$r_3$','interpreter','latex')