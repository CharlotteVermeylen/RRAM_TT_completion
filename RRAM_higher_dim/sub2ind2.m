function ndx = sub2ind2(siz,V)
    %SUB2IND Linear index from multiple subscripts.
    %   SUB2IND is used to determine the equivalent single index
    %   corresponding to a given set of subscript values.
    %
    %   IND = SUB2IND(SIZ,I,J) returns the linear index equivalent to the
    %   row and column subscripts in the arrays I and J for a matrix of
    %   size SIZ. 
    %
    %   IND = SUB2IND(SIZ,I1,I2,...,IN) returns the linear index
    %   equivalent to the N subscripts in the arrays I1,I2,...,IN for an
    %   array of size SIZ.
    %
    %   I1,I2,...,IN must have the same size, and IND will have the same size
    %   as I1,I2,...,IN. For an array A, if IND = SUB2IND(SIZE(A),I1,...,IN)),
    %   then A(IND(k))=A(I1(k),...,IN(k)) for all k.
    %
    %   Class support for inputs I,J: 
    %      float: double, single
    %      integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
    %
    %   See also IND2SUB.

    %   Copyright 1984-2015 The MathWorks, Inc.
    %   Adapted version by Charlotte.
    
    siz = double(siz);
    lensiz = length(siz);
    if lensiz < 2
        error(message('MATLAB:sub2ind:InvalidSize'));
    end 

    ndx = V(:,1)+(V(:,2)-1).*siz(1);

    k = cumprod(siz);
    for i = 3:size(V,2)
        v = V(:,i);
        %%Input checking
        if (any(min(v(:)) < 1)) || (any(max(v(:)) > siz(i)))
            %Verify subscripts are within range
            error(message('MATLAB:sub2ind:IndexOutOfRange'));
        end
        ndx = ndx + (v-1)*k(i-1);
    end
end