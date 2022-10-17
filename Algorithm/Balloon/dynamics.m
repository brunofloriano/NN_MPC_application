function next_state = dynamics(previous_state,control_input,delta_t)

    horizontal_position_x = previous_state(1);
    vertical_position_z = previous_state(2);
    u = control_input;

    velocity_gradient_alpha = 1e-3;
    noise_std_deviation = sqrt(1500);

    noise_xi = noise_std_deviation*randn(1);
    next_x = horizontal_position_x + delta_t*(velocity_gradient_alpha*vertical_position_z + noise_xi);
    next_z = u;

    next_state = [next_x;next_z];