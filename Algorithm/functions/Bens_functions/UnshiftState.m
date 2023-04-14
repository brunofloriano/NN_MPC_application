function state = UnshiftState(shift_state,d)
    state = zeros(d,1);
    for i=1:d
        if(mod(shift_state(i,1),2)==0)
            state(i,1)=shift_state(i,1)/2;
        else
            state(i,1)=(shift_state(i,1)+1)/-2;
        end
    end
end