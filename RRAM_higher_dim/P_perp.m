function P = P_perp(X)
    
    P = eye(length(X)) - X*X';

end