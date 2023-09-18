function next_state = dynamics(previous_state,control_input,delta_t)

    x = previous_state;
    u = control_input;
    
    a = 0.3;
    d = 0.8;
    A = [0 1;-d -a];
    B = [0 0; 1 0];
    next_state = previous_state + delta_t*(A*x + B*u);