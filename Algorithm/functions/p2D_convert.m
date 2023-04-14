function D = p2D_convert(p,G)
n = length(p)^2;
D.P = zeros(n+1,1); D.P(1) = 0;

D.key = zeros(n+1,1); D.key(1) = -1;
D.v = D.P; D.f = D.P; D.u = D.P; D.w = D.P; l = 2;

x = G.x; y = G.y;
xc = round(x/G.dx); yc = round(y/G.dx);

for i = 1:length(xc)
    for j = 1:length(yc)
        state = [xc(i);yc(j)];

        D.key(l,1) = state_conversion(state, G);
        D.P(l,1) = p(i,j);
        %D.P(l,1) = ((2*pi)^G.d*dist.Pdet)^(-1/2)*exp((-1/2)*(state.*G.dx-dist.mean)'*dist.Pinv*(state.*G.dx-dist.mean));
        l=l+1;
    end
end
%prob_sum = (G.dx)^G.d.*sum(D.P); D.P = D.P/prob_sum;
D.n = length(D.P);
D = Initialize_vuw(D,G);
D = boundary_conditions(D,G);
end