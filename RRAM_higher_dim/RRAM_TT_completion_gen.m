function [X,cost,test,stats,ranks] = RRAM_TT_completion_gen(A_Omega,Omega,A_Gamma,Gamma,X0,opts)
    % Riemannian rank adaptive method (RRAM) proposed in 
    % C. Vermeylen and M. Van Barel, A Riemannian rank-adaptive method for 
    % higher-order tensor completion in the tensor-train format. Arxiv (2024)
    
    X = X0;
    stats.time = 0;
    stats.rankidx = [];
    norm_A_Gamma = norm(A_Gamma);
    norm_A_Omega = norm(A_Omega);
    tic
    test = sqrt(2*func(A_Gamma,X0,Gamma))/norm_A_Gamma;
    test_current = test;
    cost = [];
    r = X0.rank;
    ranks = [];
    n = X0.size;
    grad = struct; grad.size = n; grad.sub = Omega; grad.val = euclidgrad(A_Omega,X0,Omega);
    d = X0.order;
    
    for k = 1:opts.kmax
        Xold = X;
        [X,cost_tmp,~,~] = completion_orth(A_Omega,Omega,A_Gamma,Gamma,X0,opts);
        cost_current = cost_tmp(end);
        test_old = test_current; test_current = sqrt(2*func(A_Gamma,X,Gamma))/norm_A_Gamma;
           
        progress_test = (test_current-test_old)/test_old; % preferably smaller than 0
        if  progress_test > opts.eps_gamma
            disp(['Rel. test func. increased more than ' num2str(opts.eps_gamma) '. Returning X_{k-1}.']);
            X = Xold;
            return
        end
        
        stats.time = [stats.time, toc]; 
        ranks = [ranks; r];
        stats.rankidx = [stats.rankidx, length(cost_tmp)];          
        cost = [cost; cost_tmp];
        test = [test, test_current];
        
        if cost(end) < opts.eps_omega
            disp(['Cost function converged!' cost(end) ' < ' opts.eps_omega]);
            return
        end

        if k < opts.kmax
            Xt = rank_decrease_gen(X,opts);
            rt = Xt.rank;
            if sum(r-rt) ~= 0 
                cost_trunc = sqrt(2*func(A_Omega,Xt,Omega))/norm_A_Omega;
                test_trunc = sqrt(2*func(A_Gamma,Xt,Gamma))/norm_A_Gamma;
                if (cost_trunc < cost_current) || (test_trunc < test_current)
                    X0 = Xt; rold = r; r = rt;
                    disp(['Truncated rank from ' num2str(rold) ' to ' num2str(r) ':']);
                    disp(['Rel. cost func. after rank red.:', num2str(cost_trunc)]);
                    disp(['Rel. test func. after rank red.:', num2str(test_trunc)]);
                    continue
                end
            end
            
            opts.s = ones(size(r));
            for i=1:d-1
                opts.s(i+1) = min(opts.r_max(i)-r(i+1),opts.s_max(i));
            end
            
            if sum(opts.s(2:d)) > 0
                rold = r; grad.val = euclidgrad(A_Omega,X,Omega);
                X0 = rank_incr(X,A_Omega,Omega,norm_A_Omega,A_Gamma,Gamma,norm_A_Gamma,grad,opts); 
                r = X0.rank;
                disp(['Increased rank from ' num2str(rold) ' to ' num2str(r) ':']);
                disp(['Rel. cost func. after rank increase:', ...
                    num2str(sqrt(2*func(A_Omega,X0,Omega))/norm_A_Omega)]);
                disp(['Rel. test func. after rank increase:', ...
                    num2str(sqrt(2*func(A_Gamma,X0,Gamma))/norm_A_Gamma)]);
            else
                X0 = X;
            end
        end
    end

end
