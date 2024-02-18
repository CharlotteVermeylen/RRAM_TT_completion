function Xnew = retraction_high_dim(X,Ui,Vi,i,t)

    X = orthogonalize(X,i+1);

    Xi1 = X.U{i+1};
    
    r = X.rank;
    n = X.size;
    
    Ui = reshape(Ui,r(i),n(i),[]);
    Xi = cat(3,X{i},Ui);

    if i==X.order-1
        Xi1 = [Xi1; t*Vi];
    else
        Vi = reshape(Vi,[],n(i+1),r(i+2));
        Xi1 = cat(1,Xi1,t*Vi);
    end
    
    Xnew = X;
    Xnew.U{i} = Xi;
    Xnew.U{i+1} = Xi1;
    
end