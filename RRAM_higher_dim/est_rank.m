function s = est_rank(Y,s_max)
    % Estimated rank defined in Definition 3.4 in the manuscript.
    
    sv = svd(Y); 
    
    if s_max == length(sv)
        sv = [sv;0];
    end
    
    sv = sv(1:s_max+1);
    relgap = (sv(1:end-1)-sv(2:end))./sv(1:end-1);

    s = find(relgap==max(relgap),1);

end