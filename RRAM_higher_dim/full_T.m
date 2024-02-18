function Xf = full_T(X)
    
    d = length(X);
    Xf = X{1};
    for i=2:d
       Xf = mult_T(Xf,X{i}); 
    end

end