function XL = left(X,l)

    if nargin == 1
        s = size(X);
        XL = reshape(X,s(1),[]);
    else
        XL = reshape(X,l,[]);
    end

end