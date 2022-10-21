function [x,y] = kinematics(velocity,previous_theta,radius,delta_t)

angular_velocity = velocity/radius;
next_theta = previous_theta + delta_t*angular_velocity;
[x,y] = pol2cart(next_theta,radius);

