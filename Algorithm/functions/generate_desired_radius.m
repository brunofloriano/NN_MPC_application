function desired_radius = generate_desired_radius(N, max_radius)

desired_radius{1} = 0;
radius = max_radius*rand(N);

%% Equal orbits (one agent per orbit)

for agent_counter = 2:N+1
    desired_radius{agent_counter} = radius(1,agent_counter);
end