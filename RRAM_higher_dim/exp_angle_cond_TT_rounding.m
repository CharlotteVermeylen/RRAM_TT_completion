% This script regenerates the sixth figure in the manuscript.

clear

d = 4;
ni = 10;
n = ni*ones(1,d);

nexp = 20;
angle = zeros(nexp,1);

r = [1,2,2,2,1];
rA = [1,4,4,4,1];

for i=1:nexp
    rng(i)
    A = TTeMPS_randn(rA,n);

    fprintf('exp: %i \n',i)
    
    disp('TT-rounding')
    X = truncate_TTeMPS_r(A,r);
    X_full = full(X);
    nX_full = norm(X_full(:));
    disp(norm(X_full(:)));
   
    An = A(:)/norm(A(:));
    angle(i) = An'*(X_full(:)/nX_full); 
end

%%
set(0,'defaultAxesFontSize',26)
figure; 
gca = subplot(1,2,1);
h = boxplot(angle,'Widths',0.5,'Labels',{'TT-rounding'},'Positions',1);
set(h,{'linew'},{2})
ylabel('$\Big\langle \frac{\tilde{X}}{||\tilde{X}||} , \frac{A}{||A||} \Big\rangle $','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%%
gca = subplot(1,2,2);
semilogy(angle,'b-o','linewidth',2,'Markersize',10)
omega = sqrt(prod(r(3:end))/prod(rA(3:end)));
disp(omega)
hold on
semilogy(omega*ones(1,nexp),'r--','linewidth',2,'Markersize',10)
legend('TT-rounding','$\omega$','Interpreter','Latex','location','southeast')

xlabel(strcat(sprintf('%d tensors ',nexp),' $A$'),'Interpreter','Latex')
ylabel('$\Big\langle \frac{\tilde{X}}{||\tilde{X}||} , \frac{A}{||A||} \Big\rangle $','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
axis tight
