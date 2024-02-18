function X = rank_decrease_gen(X,opts)
    % method to decrease the rank in the RRAM based on the TT-rounding
    % algorithm and the defined numerical rank in definition 6.
    
    rt = num_rank(X,opts.Delta);
    r = X.rank;
    
    if sum(r-rt) ~= 0 
        X = truncate_TTeMPS_r(X,rt);
    end
    
end
