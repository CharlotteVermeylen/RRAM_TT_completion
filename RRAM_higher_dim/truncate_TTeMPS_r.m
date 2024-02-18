function X = truncate_TTeMPS_r(X,rank)
    % This function contains our implementation of the TT-rounding
    % algorithm: Algorithm 2 in 
    % I. V. Oseledets, Tensor-train decompositions. Methods and 
    % Algorithms for Scientific Computing 33.5 (2011), pp. 2295–2317.
    % In our manuscript this algorithm is recalled as Algorithm 3.
    
    X = orthogonalize(X,1);
    
    for i = 1:X.order-1
        r = size(X.U{i});
        if length(r) == 2
            r = [r,1];
        end
        [U,S,V] = svd(reshape(X.U{i},[],r(end)));
    
        s = rank(i+1);

        U = U(:,1:s);
        V = V(:,1:s);
        S = S(1:s,1:s);
        r(end) = s;
        X.U{i} = reshape(U, r);
        X.U{i+1} = mult_T(S*V',X.U{i+1});
    end
   
    X = TTeMPS(X.U);
end