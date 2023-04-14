function state = CantorUnpair(key,d)
    state = zeros(d,1);
    for i=2:d
        z=key; w=floor(((8*z+1)^(1/2)-1)/2); t=(w^2+w)/2; 
        y=z-t; x=w-y; state(d-i+1,1)=x; state(d-i+2,1) =y;
        key=x;
    end
end