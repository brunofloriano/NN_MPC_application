%% Balloons parameters

% Choose simulation parameters
tangential_velocity = 60; % in m/s
Nr = 8; % excluding the root (mothership)
N = Nr+1; % including mothership
height = 0;
level_height = 558.3;

Agents.N = N;
Agents.Nr = Nr;

%% Time parameters
tmax = 2*24*3600/1; % in seconds
%t_change = tmax/2;
t_steps = 501 + 500;%1000;
delta_t = tmax/(t_steps-1); %1; % in seconds
t73 = tmax; % time to have a 73% decrease on certainty
alpha_time_constant = 1/t73;
time_horizon = tmax/2;

Time.tmax = tmax; 
Time.t_steps = t_steps;
Time.delta_t = delta_t;
Time.t73 = t73;
Time.alpha_time_constant = alpha_time_constant;
Time.time_horizon = time_horizon;
%% Communication parameters
comm_range = 10e3;
max_root_connections = 5;

%% Map area parameters
max_radius = 200e3;
map_steps = 61; delta_x = 2*max_radius/(map_steps - 1);
encompassed_area = zeros(map_steps);
%encompassed_area = ones(map_steps);
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
covariance = round(area_covariance_matrix/delta_x^2);
std_area = def_normal_area(map_steps,[round(map_steps/2);round(map_steps/2)],covariance);

x = -max_radius:delta_x:max_radius; y = x;
G.dx=delta_x; G.dt=0.001; G.xh=G.dx/2; G.d=2; G.dif = 1; 
G.x = x; G.y = y;

Area.map_steps = map_steps; 
Area.delta_x = delta_x; 
Area.encompassed_area = encompassed_area;
Area.encompassed_area_time = encompassed_area_time;
Area.area_radius = area_radius;
Area.total_area = total_area;
Area.total_circle_area = total_circle_area;
Area.area_covariance_matrix = area_covariance_matrix;
Area.x = x;
Area.y = y;
Area.max_radius = max_radius;
Area.max_area_pdf = max_area_pdf;
Area.covariance = covariance;
Area.std_area = std_area;

%Area.encompassed_area = def_4th_area(Area);
load('p_4th.mat')
Area.encompassed_area = D2p_convert(D);
Area.encompassed_area = Area.encompassed_area/max(max(Area.encompassed_area));
Area.initial_encompassed_area = Area.encompassed_area;
encompassed_area = Area.encompassed_area;
encompassed_area_time{1} = encompassed_area;
%% Cost function parameters
% Time
alpha_t = 1;
beta_t = 0;

% Area coverage
beta_ac = 1;

% Communication
gamma_cc = 1;
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