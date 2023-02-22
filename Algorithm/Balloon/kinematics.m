function [x,y,next_theta,radius] = kinematics(previous_theta,radius,delta_t)

velocity = radius2velocity(radius);
angular_velocity = velocity/radius;
if radius == 0
    angular_velocity = 0;
end
next_theta = previous_theta + delta_t*angular_velocity;
[x,y] = pol2cart(next_theta,radius);

