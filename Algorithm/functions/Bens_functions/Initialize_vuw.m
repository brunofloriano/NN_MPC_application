function D = Initialize_vuw(D,G)
    for l=2:D.n
        current_key = D.key(l,1); current_state = double(key_conversion(current_key,G));
           
        x=G.dx.*current_state; 
        
        x_state = x; x_state(1) = x_state(1)+G.xh; 
        r = norm(x_state); theta = atan2(x_state(2),x_state(1)); 
        v1 = -G.dif*r^3*cos(theta)*4; % Potential Bowl
        D.v(l,1) = v1;
        D.u(l,1)=min(v1, 0);
        D.w(l,1)=max(v1, 0);

        y_state = x; y_state(2) = y_state(2)+G.xh; 
        r = norm(y_state); theta = atan2(y_state(2),y_state(1)); 
        v2 = -G.dif*r^3*sin(theta)*4;
        D.v(l,2) = v2;
        D.u(l,2)=min(v2, 0);
        D.w(l,2)=max(v2, 0);
    end
end