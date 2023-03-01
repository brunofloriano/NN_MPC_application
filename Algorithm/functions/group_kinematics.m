function [individual_coord, individual_angular_position, individual_state, Ggraph] = group_kinematics(individual_angular_position, individual_state, individual_coord, time_counter, delta_t, comm_range, max_root_connections,Ggraph)

Nr = length(individual_angular_position) - 1;

for agent_counter = 1:Nr+1

    theta_position = individual_angular_position{agent_counter};
    %individual_angular_vector{agent_counter}(time_counter) = theta_position;

    radius = individual_state{agent_counter+1}(1);

    if radius < 0
        %theta_position = - theta_position;
        radius = - radius;
        individual_state{agent_counter+1}(1) = - individual_state{agent_counter+1}(1);
        %disp('Negative radius warning!')
    end
    %radius_vector{agent_counter}(time_counter) = radius;

    [x,y,theta_position,radius] = kinematics(theta_position,radius,delta_t);
    x = real(x);
    y = real(y);
    position(:,agent_counter) = [x;y];
    %[theta_position,radius] = cart2pol(x,y);

    individual_coord{agent_counter}(:,time_counter) = [x;y];
    individual_angular_position{agent_counter} = theta_position;

end
Ggraph{time_counter} = routing_protocol(position,comm_range,max_root_connections);
