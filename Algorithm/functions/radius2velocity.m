function velocity = radius2velocity(radius)

radius_max = 20e3; % in m (radius of max velocity
velocity_max = 60; % in m/s
bs = 1.8;
x = 0.5;

% if radius < radius_max
%     velocity = velocity_max*(radius/radius_max);
% else
%     velocity = velocity_max*(radius_max/radius)^(0.5);
% end

velocity = velocity_max*( (radius_max/radius)^bs * exp( 1 - (radius_max/radius)^bs ) )^x;