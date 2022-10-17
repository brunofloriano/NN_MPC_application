function [new_individual_state,new_L] = add_new_agent(individual_state,new_agent_state,L,connection_index)
    N_old = lenght(individual_state);
    N_new = N_old + 1;
    new_agent_index = N_new;
    
    new_individual_state = individual_state;
    new_individual_state{N_new} = new_agent_state;
    
    new_L = [L zeros(N_old,1);zeros(1,N_old) 0];
    
    new_L(conection_index,new_agent_index) = -1;
    new_L(new_agent_index,conection_index) = -1;
    
%     for agents_counter = 1:N_old
%         distance = norm(individual_state{agents_counter} - new_agent_state);
%         if sum(agents_counter == connection_index)  %distance < 5
%             laplacian_component = -1;
%         else
%             laplacian_component = 0;
%         end
%         new_L(agents_counter,new_agent_index) = laplacian_component;
%         new_L(new_agent_index,agents_counter) = laplacian_component;
%     end
    
    new_L(new_agent_index,new_agent_index) = -sum(new_L(new_agent_index,:));
    
end