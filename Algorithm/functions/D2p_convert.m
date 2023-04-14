function p = D2p_convert(D)
n = D.n;
n2 = round(sqrt(n-1));
p = zeros(n2);
l = 2;

for i = 1:n2
    for j = 1:n2
        p(i,j) = D.P(l,1);
       l=l+1;
    end
end

end