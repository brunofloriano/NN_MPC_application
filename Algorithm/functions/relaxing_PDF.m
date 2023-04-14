function D = relaxing_PDF(D,G)

%D = p2D_convert(p,G);

%K = RHS_P(D,G);

k1=RHS_P(D,G); D.P = D.P + (G.dt/2).*k1;
k2=RHS_P(D,G); D.P = D.P + (G.dt/2).*k2;
k3=RHS_P(D,G); D.P = D.P + G.dt.*k3;
k4=RHS_P(D,G);
D_new.P=D.P+(G.dt/6).*k1+(G.dt/3).*(k2+k3)+(G.dt/6).*k4;
diff = sum(abs(D.P - D_new.P));
D.P = D_new.P;

%p = D2p_convert(D);