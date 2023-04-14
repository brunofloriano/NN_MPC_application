function state = key_conversion(key,G)
    if(G.d==1)
        shift_state = [key]; 
    else
        shift_state = CantorUnpair(key,G.d);
    end
    state = UnshiftState(shift_state,G.d);
end