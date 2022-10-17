function individual_state = group_dynamics(individual_state,K,Adjacency,delta_t)
    N = length(individual_state) - 1;
    m = size(K,1);
    for agent_counter_i = 1:N
        weighted_sum = 0;
        for agent_counter_j = 1:N
            epsilon_ij{agent_counter_j} = individual_state{agent_counter_j+1} - individual_state{agent_counter_i+1};
            weighted_sum = weighted_sum + Adjacency(agent_counter_i,agent_counter_j)*epsilon_ij{agent_counter_j};
        end

        individual_control_input = -K*weighted_sum;
        individual_state{agent_counter_i+1} = dynamics(individual_state{agent_counter_i+1},individual_control_input,delta_t);
    end
    individual_state{1} = dynamics(individual_state{1},zeros(m,1),delta_t); % Virtual agent dynamics
end