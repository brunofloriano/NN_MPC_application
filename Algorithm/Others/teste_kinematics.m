clear all; close all; clc;

delta_t = 0.1;
tmax = 1e5;
t = 0:delta_t:tmax;
v = 80/3.6;
theta = zeros(length(t),1);
R = 40e3;

omega = v/R;

for time_counter = 2:length(t)
    theta(time_counter) = theta(time_counter-1) + delta_t*omega;
end

[x,y] = pol2cart(theta,R*ones(length(theta),1));

plot(x,y)