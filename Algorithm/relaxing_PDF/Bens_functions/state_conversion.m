function key = state_conversion(state,G)
    shift_state = ShiftState(state,G.d);
    if (G.d==1)
        key = shift_state(1);
    else
        key = CantorPair(shift_state);
    end
end