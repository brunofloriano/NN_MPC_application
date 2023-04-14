function key = CantorPair(state)
    if(length(state)>2)
        last = state(end,1); state(end,1) = [];
        x = CantorPair(state); y = last;
        key = (1/2)*(x+y)*(x+y+1)+y;
    else
        x=state(1,1); y=state(2,1);
        key = (1/2)*(x+y)*(x+y+1)+y;
    end
end