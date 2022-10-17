function next_state = dynamics(previous_state,control_input,delta_t)
    
    L = 0.3;
    
    x = previous_state;
    xf = x(1);
    yf = x(2);
    theta = x(3);
    delta = x(4);
    
    u = control_input(1);
    vf = control_input(2);
    
    f1 = vf*cos(theta+delta);
    f2 = vf*sin(theta+delta);
    f3 = vf*sin(delta)/L;
    f4 = u;
    
    f = [f1 f2 f3 f4]';

    next_state = previous_state + delta_t*f;