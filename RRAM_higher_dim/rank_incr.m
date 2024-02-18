function X = rank_incr(X,A_Omega,Omega,norm_A_Omega,A_Gamma,Gamma,norm_A_Gamma,grad,opts)
    % Method to increase the rank in the RRAM: Algorithm 2 in the
    % manuscript.
    d = X.order;
    r = X.rank;
    n = X.size;
    f_Omega = sqrt(2*func(A_Omega,X,Omega))/norm_A_Omega;
    f_Gamma = sqrt(2*func(A_Gamma,X,Gamma))/norm_A_Gamma;

    for i=2:d
        grad.val = euclidgrad(A_Omega,X,Omega); 
        Pi = PiXY(X,grad,i-1); 
        smax = min([opts.s(i),rank(Pi,10^-13)]);
        si = est_rank(Pi,smax);
        rnewi = r(i) + si;
        [Ui,Si,Vi] = svd(Pi);
        Ui = Ui(:,1:si); Si = Si(1:si,1:si); Vi = Vi(:,1:si);
        Vi = Si*Vi';
        PX = X;
        PX.U{i-1} = reshape(Ui,r(i-1),n(i-1),si);
        PX.U{i} = reshape(Vi,si,n(i),r(i+1));
        
        dx_Omega = PX(grad.sub);
        t = -dx_Omega'*grad.val/norm(dx_Omega)^2;
        fprintf('t: %g',t)
        Xnew = retraction_high_dim(X,Ui,Vi,i-1,t);
        
        f_Omega_new = sqrt(2*func(A_Omega,Xnew,Omega))/norm_A_Omega;
        f_Gamma_new = sqrt(2*func(A_Gamma,Xnew,Gamma))/norm_A_Gamma;
        disp(['Trying to increase rank r(' num2str(i) ') from ' num2str(r(i))...
            ' to ' num2str(rnewi) ':']);
        disp(['Rel. cost func. after rank increase:', num2str(f_Omega_new)]);
        disp(['Rel. test func. after rank increase:', num2str(f_Gamma_new)]);
        disp(f_Omega-f_Omega_new); 
        disp(f_Gamma-f_Gamma_new);
        if (f_Omega-f_Omega_new) > 10^-10 && (f_Gamma-f_Gamma_new) > 10^-10
            X = Xnew;
            r = X.rank;
            f_Omega = f_Omega_new;
            f_Gamma = f_Gamma_new;
        end
    end  
end