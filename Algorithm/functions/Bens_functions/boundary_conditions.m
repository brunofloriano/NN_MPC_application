function D = boundary_conditions(D,G)
    for l=2:D.n
        key = D.key(l,1); state = key_conversion(key,G); pos = G.dx.*state; 
        if(max(abs(pos))==3)
            D.P(l,1) = 0;
        end
    end     
end