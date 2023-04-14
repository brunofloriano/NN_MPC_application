function K = RHS_P(D,G)
    diffusion = zeros(D.n,1); K = zeros(D.n,1);

    for l=2:D.n %Flux out
        l_key = D.key(l,1); l_state = key_conversion(l_key,G); 
        for q=1:G.d
            k_state = l_state; k_state(q) = k_state(q)+1; k_key = state_conversion(k_state,G);
            i_state = l_state; i_state(q) = i_state(q)-1; i_key = state_conversion(i_state,G);
            k_pos = find(D.key==k_key); i_pos = find(D.key==i_key); 
            if(isempty(k_pos)),k_pos = 1;end
            if(isempty(i_pos)),i_pos = 1;end
            D.f(l,q) = D.u(l,q)*D.P(k_pos,1) + D.w(l,q)*D.P(l,1);
            diffusion(l,1) = diffusion(l,1) + D.P(k_pos,1)-2*D.P(l,1)+D.P(i_pos,1); 
        end
    end
    
    for d=1:G.d %Flux out w/ CTU
        for l=2:D.n
            l_key = D.key(l,1); l_state = key_conversion(l_key,G); 
            i_state = l_state; i_state(q) = i_state(q)-1; i_key = state_conversion(i_state,G);
            i_pos = find(D.key==i_key); if(isempty(i_pos)),i_pos = 1;end

            F = G.dt*(D.P(l,1)-D.P(i_pos,1))/(2*G.dx);
            for e=1:G.d
                if e~=d
                    j_state = l_state; j_state(e) = j_state(e)-1; j_key = state_conversion(j_state,G);
                    j_pos = find(D.key==j_key); if(isempty(j_pos)),j_pos = 1;end
                    p_state = i_state; p_state(e) = p_state(e)-1; p_key = state_conversion(p_state,G);
                    p_pos = find(D.key==p_key); if(isempty(p_pos)),p_pos = 1;end

                    D.f(l,e) = D.f(l,e)-D.w(l,e)*D.w(i_pos,d)*F;
                    D.f(j_pos,e) = D.f(j_pos,e)-D.u(j_pos,e)*D.w(i_pos,d)*F;
                    D.f(i_pos,e) = D.f(i_pos,e)-D.w(i_pos,e)*D.u(i_pos,d)*F;
                    D.f(p_pos,e) = D.f(p_pos,e)-D.u(p_pos,e)*D.u(i_pos,d)*F;
                end
            end

            if(D.v(i_pos,d) > 0)
                i_i_state = i_state; i_i_state(d) = i_i_state(d)-1; i_i_key = state_conversion(i_i_state,G);
                i_i_pos = find(D.key==i_i_key); if(isempty(i_i_pos)),i_i_pos = 1;end
                th = (D.P(i_pos,1) - D.P(i_i_pos,1))/(D.P(l,1)-D.P(i_pos,1));
            else
                k_state = l_state; k_state(q) = k_state(q)+1; k_key = state_conversion(k_state,G);
                k_pos = find(D.key==k_key); if(isempty(k_pos)),k_pos = 1;end
                th = (D.P(k_pos,1) - D.P(l,1))/(D.P(l,1)-D.P(i_pos,1));
            end

            t = abs(D.v(i_pos,d));
            D.f(i_pos,d) = D.f(i_pos,d) + t*(G.dx/G.dt-t)*F*MC(th);
        end
    end
          
    for l=2:D.n
        l_key = D.key(l,1); l_state = key_conversion(l_key,G); l_pos = abs(G.dx.*l_state);
        if (max(l_pos)~=3)
            K(l,1) = K(l,1) + (1/G.dx^2)*(G.dif*diffusion(l,1));
            for d=1:G.d
                i_state = l_state; i_state(d) = i_state(d)-1; i_key = state_conversion(i_state,G);
                i_pos = find(D.key==i_key); 
                if(isempty(i_pos))
                    K(l,1) = K(l,1)-(1/G.dx)*(D.f(l,d))+(1/G.dx^2)*(G.dif*diffusion(l,1));
                else
                    K(l,1) = K(l,1)-(1/G.dx)*(D.f(l,d)-D.f(i_pos,d));
                end
            end
        end
    end
end