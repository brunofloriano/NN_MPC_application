h = 0.8394; %LMI minimization

desired_individual_state{1} = [0 0]';
desired_individual_state{2} = [0 0]';
desired_individual_state{3} = [0 0]';
desired_individual_state{4} = [0 0]';

desired_individual_state{5} = [0 0]'; % agent to be added later

% WW20 = 1 => WW20
% WW20 = 0 => SS16
ww20 = 0;

nNeurons = 4;
delta_t = 0.1;

% INDIVIDUAL SYSTEM DYNAMICS
a = 0.3;
d = 0.8;
c = d;
alpha = [1 2 3]; 

if ww20 == 1
    phi = [1 3 6]; %WW20
    phi_dot = [4 3 -2]; %WW20
    A = [0 1;-d -a]; %WW20
    tmax = 30;
    delta = 0.3;
else
    phi = [0 4 7]; %SS16
    phi_dot = [3 2 -1]; %SS16
    A = [0 1;0 -1]; %SS16
    tmax = 50; %18;
    delta = 0.2;
    K = [1 0;0 0];
end

B = [1];
t_change = 10;
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
L{1} = [1 0 -1 0;0 0 0 0;0 0 0 0;0 0 0 0];
L{2} = [0 0 0 0;0 1 -1 0;0 -1 1 0;0 0 0 0];

L{3} = [1 0 0 -1;0 1 -1 0;0 -1 2 -1;-1 0 -1 2];
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