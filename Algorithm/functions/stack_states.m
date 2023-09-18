function global_state = stack_states(individual_state)
N = length(individual_state) - 1;
n = length(individual_state{1});

global_state = zeros(n*N,1);
for cnt = 1:N
    vector_start = (cnt-1)*n + 1;
    vector_end = cnt*n;
    global_state(vector_start:vector_end,1) = individual_state{cnt+1};
end