function Y = struct2full(Ys)

    Y = zeros(Ys.size);
    Y(sub2ind2(Ys.size,Ys.sub)) = Ys.val;

end