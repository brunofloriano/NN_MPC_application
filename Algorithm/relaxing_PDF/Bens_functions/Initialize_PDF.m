function D = Initialize_PDF(G,dist)
    D.P = zeros(dist.n+1,1); D.P(1) = 0; D.key = zeros(dist.n+1,1); D.key(1) = -1;
    D.v = D.P; D.f = D.P; D.u = D.P; D.w = D.P; l = 2;  

    for i=round((dist.mean(1,1)-3*dist.std(1,1))/G.dx):round((dist.mean(1,1)+3*dist.std(1,1))/G.dx)
        for j=round((dist.mean(2,1)-3*dist.std(2,1))/G.dx):round((dist.mean(2,1)+3*dist.std(2,1))/G.dx)
            state = [i;j]; D.key(l,1) = state_conversion(state, G);  
            D.P(l,1) = ((2*pi)^G.d*dist.Pdet)^(-1/2)*exp((-1/2)*(state.*G.dx-dist.mean)'*dist.Pinv*(state.*G.dx-dist.mean));
            l=l+1;
        end
    end
    prob_sum = (G.dx)^G.d.*sum(D.P); D.P = D.P/prob_sum; D.n = length(D.P); 
end