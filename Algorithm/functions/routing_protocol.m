function Ggraph = routing_protocol(position,comm_range,max_root_connections)

N = length(position)-1;
root_position = position(:,1);

Ggraph = graph([],[],[],N+1);

root_connections = 0;
dist_root = zeros(1,N);
for agent_counter_i = 1:N
    dist_root(agent_counter_i) = norm( position(:,agent_counter_i+1) -  root_position);
end

[dist_root_sorted,sorted_index] = sort(dist_root);

for agent_counter_i = sorted_index %2:N+1
    dist_root_i = dist_root(agent_counter_i);

    if dist_root_i < comm_range && root_connections < max_root_connections
        Ggraph = addedge(Ggraph,1,agent_counter_i+1);
        %plot(Ggraph,XData=position(1,:),YData=position(2,:))
        root_connections = root_connections + 1;
    else

        %neighbors_set = [];
        min_cost = 1e100;
        min_cost_index = -1;
        for agent_counter_j = sorted_index %1:N
            dist_ij = norm( position(:,agent_counter_i+1) -  position(:,agent_counter_j+1) );
            dist_root_j = dist_root(agent_counter_j);

            if agent_counter_j ~= agent_counter_i && dist_ij < comm_range
                %neighbors_set = [neighbors_set agent_counter_j];
                cost = dist_ij;
                if dist_root_j < dist_root_i && cost < min_cost
                    min_cost = cost;
                    min_cost_index = agent_counter_j;
                end
            end
        end

        if min_cost_index > 0
            Ggraph = addedge(Ggraph,min_cost_index+1,agent_counter_i+1);
            %plot(Ggraph,XData=position(1,:),YData=position(2,:))
        end
    end
end