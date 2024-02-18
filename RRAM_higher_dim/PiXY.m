function Pi = PiXY(X,Y,i)
    % This function computes the matrices $P_i(X,Y)$ defined in the manuscript
    % C. Vermeylen and M. Van Barel, A Riemannian rank-adaptive method for 
    % higher-order tensor completion in the tensor-train format, ArXiv
    % (2024),
    % in equation (32). 
    
    n = X.size;
    d = length(n);
    r = X.rank;
    X1 = orthogonalize(X,1);
    Xd = orthogonalize(X,d);
    
    if isstruct(Y)
        Y = struct2full(Y);
    end
    Yi = reshape(Y,prod(n(1:i-1)),n(i)*n(i+1),prod(n(i+2:end)));
    
    P_Xi_perp = eye(n(i)*r(i)) - right(Xd.U{i},r(i+1))*right(Xd.U{i},r(i+1))'; % maybe not necessary
    P_Xi1_perp = eye(n(i+1)*r(i+2)) - left(X1.U{i+1},r(i+1))'*left(X1.U{i+1},r(i+1)); % maybe not necessary
    
    if i==1
        XLi = left(full_T(X1.U(i+2:end)));
        Pi = mult_T(Yi,XLi');
    elseif i == d-1
        XRi = right(full_T(Xd.U(1:i-1)),r(i));
        Pi = mult_T(XRi',Yi);
    else
        XLi = left(full_T(X1.U(i+2:end)));
        XRi = right(full_T(Xd.U(1:i-1)),r(i));
        Pi = mult_T(mult_T(XRi',Yi),XLi');
    end
    
    Pi = P_Xi_perp*reshape(Pi,r(i)*n(i),[])*P_Xi1_perp;
end