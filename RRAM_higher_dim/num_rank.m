function r = num_rank(X,Delta)
    % Numerical rank defined in the manuscript in definition 6.
    
    X = orthogonalize(X,1);

    ranks = X.rank;
    r = ones(size(ranks));
    
    for i = 1:X.order-1
        [~,S,V] = svd(right(X.U{i},ranks(i+1)),'econ');
        
        sv = diag(S);
        if length(sv) == 1
            s = 1;
        else
            rel_gap = (sv(1:end-1)-sv(2:end))./sv(1:end-1);
            I = find(rel_gap > Delta);
            if isempty(I)
                s = ranks(i+1);
            else
                s = I(1);
            end
        end
        
        r(i+1) = s;
        X.U{i+1} = mult_T(S*V',X.U{i+1});
    end

end
