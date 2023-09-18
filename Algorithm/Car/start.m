% INDIVIDUAL SYSTEM DYNAMICS
alpha = [1 2 3];
desired_individual_state{1} = [0 0 0 0]';
desired_individual_state{2} = [-17 16 0 0]';
desired_individual_state{3} = [9 -18 0 0]';
desired_individual_state{4} = [25 -2 0 0]';

delta_t = 0.1;

phi = 5*[-0.6 1 2.5 -1];
phi_dot = 5*[0.6 -1 0.5 1];
theta_dot = pi/2*[1 0 -1 0];
delta_dot = [0 0 0 0];
A = [-1 -1;0 0];
B = [1 1;1 1];
% Bw = [1 0;0 1];
tmax = 60;
nNeurons = 20;

% COMMUNICATION GRAPH
L{1} = [1 0 -1;0 0 0;0 0 0];
L{2} = [0 0 0;0 1 -1;0 -1 1];

Delta = 0.01;
delta = 0;


Pi_estimated = [0.95    0.05;    0.02    0.98]; %Pi_1
%Pi_estimated =[0.6 0.4;0.3 0.7]; %Pi_2
%Pi_estimated =[0.5 0.5;0.5 0.5]; %Pi_3
%Pi_estimated =[0.05 0.95;0.98 0.02]; %Pi_4
mu = 0;
tau = 0;

% CONTROL SYSTEM
%K = [0.4683 -0.2158;-0.3932 0.3281];
%K = [1 0;0 0]; %rng(300); % For reproducibility
Km{1} = [0.8955 1.719];
Km{2} = [1.7479 3.2955];

% SAVE NAME
save_mainname = 'dataCarnKMPC ';