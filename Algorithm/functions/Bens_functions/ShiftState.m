function shift_state = ShiftState(state,d)
    shift_state = zeros(d,1);
    for i=1:d
        if(state(i,1)<0)
            shift_state(i,1)=-2*state(i,1)-1;
        else
            shift_state(i,1)=2*state(i,1);
        end
    end
end