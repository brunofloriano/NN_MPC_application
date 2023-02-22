% function individual_coord = Emax(individual_coord,individual_angular_position,individual_state,delta_t,time_counter)
% 
% N = length(individual_coord);
% 
% for agent_counter = 1:N
% 
%     theta = individual_angular_position{agent_counter};
%     radius = individual_state{agent_counter}(1);
% 
%     [x,y] = kinematics(theta,radius,delta_t);
%     position(:,agent_counter) = [real(x);real(y)];
%     [theta,radius] = cart2pol(real(x),real(y));
% 
%     individual_coord{agent_counter}(:,time_counter) = [real(x);real(y)];
%     individual_angular_position{agent_counter} = theta;
% 
% end

if strcmp(save_mainname,'dataBalloon ')

    for agent_counter = 1:Nr+1

        theta_position = individual_angular_position{agent_counter};
        %individual_angular_vector{agent_counter}(time_counter) = theta_position;

        radius = individual_state{agent_counter+1}(1);
        
        if radius < 0
            %theta_position = - theta_position;
            radius = - radius;
            individual_state{agent_counter+1}(1) = - individual_state{agent_counter+1}(1);
            disp('Negative radius warning!')
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

end