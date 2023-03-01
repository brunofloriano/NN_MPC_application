%% Balloons parameters

% Choose simulation parameters
tmax = 24*3600; % in seconds
%t_change = tmax/2;
delta_t = tmax/1000; %1; % in seconds
t73 = tmax; % time to have a 73% decrease on certainty
alpha_time_constant = 1/t73;
time_horizon = tmax/2;

tangential_velocity = 60; % in m/s
max_radius = 200e3;
Nr = 8; % excluding the root (mothership)
N = Nr+1; % including mothership
height = 0;
level_height = 558.3;

%% Communication parameters
comm_range = 10e3;
max_root_connections = 5;

%% Map parameters
map_steps = 200;
encompassed_area = zeros(map_steps);
encompassed_area_time{1} = encompassed_area;
area_radius = 10e3;
total_area = pi*max_radius^2;
total_circle_area = def_circle_area(map_steps);
area_normal_variance = area_radius^2;
area_normal_covariance = 0; %area_normal_variance/2;
area_covariance_matrix = area_normal_variance*eye(2) + area_normal_covariance*[0 1;1 0];
%plane_size = max_radius;
%plane_step = 2*plane_size/map_steps;
%covariance = round(area_covariance_matrix/plane_step^2);
%normal_area = def_normal_area(map_steps,[round(map_steps/2);round(map_steps/2)],covariance);
max_area_pdf = (2*pi)^(-1)*(det(area_covariance_matrix))^(-1/2);

%% Cost function parameters
alpha_t = 1;
beta_t = 0;

beta_ac = 1;

gamma_cc = 1    ;
beta_cc = 1/(Nr + gamma_cc*Nr^2); %1;
%% NN MPC Paraemeters
h = 0.8394; %LMI minimization

for agent_counter = 1:N
    desired_individual_state{agent_counter} = [0 0]';
end

nNeurons = 6;

% INDIVIDUAL SYSTEM DYNAMICS
A = zeros(2);
B = ones(1,N);


% COMMUNICATION GRAPH
% L{1} = [1 0 -1;0 0 0;0 0 0];
% L{2} = [0 0 0;0 1 -1;0 -1 1];
% Delta = 0.01;
% 
% %Pi_estimated =[-1 1;1 -1]; %[0.9905    0.0095;    0.0097    0.9903];
% %Pi_estimated =[-5 5;2 -2];
% %Pi_estimated =[0.6 0.4;0.3 0.7];
% %Pi_estimated =[0.5 0.5;0.5 0.5];
% Pi_estimated =[0.05 0.95;0.98 0.02];

% ADD AND REMOVE AGENTS
L{1} = zeros(N);
L{2} = L{1};
L{3} = L{2};
L_cumulative = L{1};
%L{1} = [1 0 -1 0;0 0 0 0;0 0 0 0;0 0 0 0];
%L{2} = [0 0 0 0;0 1 -1 0;0 -1 1 0;0 0 0 0];
%L{3} = [1 0 0 -1;0 1 -1 0;0 -1 2 -1;-1 0 -1 2];

Delta = 0.01;

%Pi_estimated =[-1 1;1 -1]; %[0.9905    0.0095;    0.0097    0.9903];
%Pi_estimated =[-5 5;2 -2];
%Pi_estimated =[0.6 0.4;0.3 0.7];
%Pi_estimated =[0.5 0.5;0.5 0.5];
Pi_estimated =[0.05 0.95 0;0.98 0.02 0;1 0 0];

mu = 0;
tau = 0;

% CONTROL SYSTEM
%K = [0.4683 -0.2158;-0.3932 0.3281];
%K = [1 0;0 0]; %rng(300); % For reproducibility
%K = [0 0;0 0];

% SAVE NAME
save_mainname = 'dataBalloon ';