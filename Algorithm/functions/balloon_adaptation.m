if strcmp(save_mainname,'dataBalloon ')
    Ggraph{1} = routing_protocol(position,comm_range,max_root_connections);

    L{1} = laplacian(Ggraph{1});
    L{3} = -1*ones(N) + (N)*eye(N); % Complete graph

    nr = 1;
    desired_individual_state{1} = individual_state{1};

    %initial_state_augmented = [];
    teste = initial_state;
    for agent_counter = 1:Nr+1
        initial_state{agent_counter+1} = teste{agent_counter};
        %initial_state_augmented = [initial_state_augmented;teste{agent_counter}];
    end
    desired_individual_state = initial_state;
end