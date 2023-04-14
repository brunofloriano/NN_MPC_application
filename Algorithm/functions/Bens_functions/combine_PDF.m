function D_p = combine_PDF(D_p, D_ag,G)
    for l=2:D_p.n
        current_key = D_p.key(l,1); ag_pos = find(D_ag.key==current_key); 
        if(~isempty(ag_pos))
            D_p.P(l,1) = max(D_p.P(l,1) - D_ag.P(ag_pos,1),0);
        end
    end 
    prob_sum = (G.dx)^G.d.*sum(D_p.P); D_p.P(3:D_p.n-1) = D_p.P(3:D_p.n-1)/prob_sum;
end