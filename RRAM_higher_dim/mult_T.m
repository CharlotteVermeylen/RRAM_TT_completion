function C = mult_T(A,B)
    % Computes the tensor product A*B = A^R*B^L.
    
    sa = size(A);
    sb = size(B);
    
    if sb(1) == 1 && sa(end) ~= 1
        sa = [sa, 1];
    end

    C = reshape(A,[],sa(end))*left(B);

    sa(end) = [];
    sb(1) = [];
    s = [sa,sb];
    
    C = reshape(C,s);
end