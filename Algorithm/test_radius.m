addpath functions

r = 0:0.1e3:500e3;
v = 0*r;

for counter = 1:length(r)
    v(counter) = radius2velocity(r(counter));
end

plot(r/1e3,v)