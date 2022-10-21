% This code implements a group of balloons around a hurricane
% by Bruno R. O. Floriano

clear; close all; clc;
addpath 'functions'

% Include the system functions and parameters
addpath 'Balloon'
start
save_folder = 'results/';

% Choose simulation parameters
delta_t = 0.1; % in seconds
tmax = 1*3600; % in seconds
tangential_velocity = 60; % in m/s

SAVE_DATA = 0;

% Generate initial conditions
N = length(initial_conditions);

initial_angular_position = initial_angular_positions();
initial_condition = initial_conditions();

individual_angular_position = initial_angular_position;
individual_state = initial_condition;

for agent_counter = 1:N
    [x,y] = pol2cart(individual_angular_position{agent_counter},individual_state{agent_counter}(1));
    individual_coord{agent_counter}(:,1) = [x;y];
end

% Simulation
t = 0:delta_t:tmax-delta_t;

% Time loop
for time_counter = 2:length(t)
    for agent_counter = 1:N
        individual_state{agent_counter} = dynamics(individual_state{agent_counter},0,delta_t);
        theta = individual_angular_position{agent_counter};
        radius = individual_state{agent_counter}(1);

        [x,y] = kinematics(tangential_velocity,theta,radius,delta_t);
        [theta,radius] = cart2pol(x,y);

        individual_coord{agent_counter}(:,time_counter) = [x;y];
        individual_angular_position{agent_counter} = theta;
    end
end

for agent_counter = 1:N
    plot(individual_coord{agent_counter}(1,:),individual_coord{agent_counter}(2,:))
    hold on
end



% Save data
if SAVE_DATA == 1
    comment = ['This simulation was perform with:' newline];
    comment = [comment '']; %Include here what is different for this sim
    clocktime = clock;
    time_now = [num2str(clocktime(1)) '-' num2str(clocktime(2)) '-' num2str(clocktime(3)) '-' num2str(clocktime(4)) '-' num2str(clocktime(5))];
    savefolderlocal = [save_folder '/' save_mainname time_now];
    save_file = [savefolderlocal '/' save_mainname  time_now '.mat'];

    mkdir(savefolderlocal);
    save(save_file);
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '.eps']);
    saveas(gcf,[savefolderlocal '/' save_mainname  time_now '.fig']);
end