function energy_increment = time_energy_increment(fk, alpha_time_constant, t, time_horizon,delta_t)

time_cycle = 0; %floor(time/time_horizon);
shifted_time = t - time_horizon*time_cycle;

energy_increment = fk - alpha_time_constant*exp(-alpha_time_constant*shifted_time)*delta_t;