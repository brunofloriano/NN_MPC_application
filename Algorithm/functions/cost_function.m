function J = cost_function(x,L)
    P = 0.5*(L+L');
    J = x'*P*x;
end