function XR = right(X,r2)

 if nargin == 1
     r = size(X);
     XR = reshape(X,[],r(end));
 else 
     XR = reshape(X,[],r2);
 end

end