clear

d = 5;
nn = 10;
n = nn*ones(1,d);

rA = [1,3*ones(1,d-1),1];
rng(1)
A = TTeMPS_randn(rA, n);
A_full = full(A); 

rx = [1, 1*ones(1,d-1), 1];
X0 = TTeMPS_randn(rx, n);

sizeOmega = round(0.3*prod(n));

rng(1)
Omega_ind = randperm(prod(n),sizeOmega)';
Omega = ind2sub2(n,Omega_ind);

sizeOmega_C = []; 
Omega_C = cell(1,d);
for i=1:d
   Omega_C{i} = [];
end
Omega_C_ind = []; A_Omega_C = [];

if isnumeric(A)
    A_Omega = A(Omega_ind);
else
    A_Omega = A(Omega);
end

%%
opts = struct('maxiter',15,'tol',10^-16,'gradtol',10^(-11),'cg',1,'verbose',0);
opts.reltol = 0;

[X,cost_tmp,test,stats] = completion_orth(A_Omega, Omega, A_Omega_C, Omega_C, X0, opts);

%%
X = orthogonalize(X,2);

A_Omega_full = zeros(n);
A_Omega_full(Omega_ind) = A_Omega;

grad = struct; grad.size = n; grad.sub = Omega; % Charlotte
grad.ind = Omega_ind;
grad.val = euclidgrad(A_Omega,X,Omega);

P = cell(1,d-1); 
for i=1:d-1
    P{i} = PiXY(X,grad,i);
end

%%
f = 21;
figure
smax1 = max(rA) +4;
tiledlayout(d-1,5, 'Padding', 'none', 'TileSpacing', 'compact');
set(0,'defaultAxesFontSize',21)
for i=1:d-1
    nexttile    
    s = svd(P{i}); s = s(1:smax1);
    p = plot(s,'-o','linewidth',2,'Markersize',10);
    ax = gca;
    ax.FontSize = f; 
    axis tight
    xlabel('$i$','Interpreter','latex','FontSize',21)
    ylabel('$\sigma_i$','Interpreter','latex','FontSize',21)
    title(strcat(sprintf('$P_{%d}',i),'(X^*,\nabla f_{\Omega}(X^*))$'),...
        'Interpreter','latex','FontSize',23)
    x = xlim(gca);
    y = ylim(gca);
    xticks(x(1):x(2));
    yticks((10.*(ceil(y(1)/10):4:floor(y(2)/10)))+10);
    
    nexttile  
    relgap = (s(1:end-1)-s(2:end))./s(1:end-1);
    plot(relgap,'-o','linewidth',2,'Markersize',10)
    axis tight
    ax = gca;
    ax.FontSize = 21; 
    xlabel('$i$','Interpreter','latex','FontSize',21)
    ylabel('$\frac{\sigma_i-\sigma_{i+1}}{\sigma_i}$','Interpreter',...
        'latex','FontSize',23)
    x = xlim(gca);
    y = ylim(gca);
    xticks(x(1):x(2));
    yticks((0.1.*(ceil(y(1)*10):2:floor(y(2)*10))));

    nexttile    
    s = svd(reshape(A_Omega_full,prod(n(1:i)),[])); s = s(1:smax1);
    plot(s(1:smax1),'-o','linewidth',2,'Markersize',10)
    axis tight
    ax = gca;
    ax.FontSize = 21; 
    xlabel('$i$','Interpreter','latex','FontSize',21)
    ylabel('$\sigma_i$','Interpreter','latex','FontSize',21)
    title(strcat('$A_\Omega',sprintf('^{<%d>}$',i)),...
        'Interpreter','latex','FontSize',23)
    x = xlim(gca);
    y = ylim(gca);
    xticks(x(1):x(2));
    yticks((100.*(ceil(y(1)/100):1:floor(y(2)/100))));
    
    nexttile  
    relgap = (s(1:end-1)-s(2:end))./s(1:end-1);
    plot(relgap,'-o','linewidth',2,'Markersize',10)
    axis tight
    ax = gca;
    ax.FontSize = 21; 
    xlabel('$i$','Interpreter','latex','FontSize',21)
    ylabel('$\frac{\sigma_i-\sigma_{i+1}}{\sigma_i}$','Interpreter',...
        'latex','FontSize',23)
    x = xlim(gca);
    y = ylim(gca);
    xticks(x(1):x(2));
    yticks((0.1.*(ceil(y(1)*10):1:floor(y(2)*10))));
    
    nexttile    
    s = svd(reshape(A_full,prod(n(1:i)),[]));
    semilogy(s(1:smax1),'-o','linewidth',2,'Markersize',10)
    axis tight
    ax = gca;
    ax.FontSize = 21; 
    xlabel('$i$','Interpreter','latex','FontSize',21)
    ylabel('$\sigma_i$','Interpreter','latex','FontSize',21)
    title(strcat('$A',sprintf('^{<%d>}$',i)),...
        'Interpreter','latex','FontSize',23)
    y = ylim(gca);
    yticks(10.^(ceil(log10(y(1))):5:floor(log10(y(2)))));
    x = xlim(gca);
    xticks(x(1):x(2));
end
