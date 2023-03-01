function energy = time_energy(time_constant, time, time_horizon)

time_cycle = 0; %floor(time/time_horizon);

shifted_time = time - time_horizon*time_cycle;

energy = exp(-time_constant*shifted_time);