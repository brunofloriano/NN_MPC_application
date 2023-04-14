function p = initial_test

x = -30:30;
y = x;

for i = 1:61

    for j = 1:61
        if norm([x(i) y(j)]) < 2
            p(i,j) = 1;
        else
            p(i,j) = 0;
        end
    end
end