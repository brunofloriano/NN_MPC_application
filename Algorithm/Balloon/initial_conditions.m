function initial_state = initial_conditions()

%    initial_state{1} = [1;4]; % Virtual agent
%    initial_state{2} = [1;4];
%    initial_state{3} = [3;3];
%    initial_state{4} = [6;-2];
%
%    initial_state{5} = [10;-7]; % agent to be added later

inner_radius = 40e3;
outter_radius = 80e3;
height = 0;

initial_state{1} = [0,height];

initial_state{2} = [inner_radius,height];
initial_state{3} = [inner_radius,height];
initial_state{4} = [inner_radius,height];
initial_state{5} = [inner_radius,height];

initial_state{6} = [outter_radius,height];
initial_state{7} = [outter_radius,height];
initial_state{8} = [outter_radius,height];
initial_state{9} = [outter_radius,height];

% initial_state{1} = [inner_radius,height];
% initial_state{2} = [inner_radius,height];
% initial_state{3} = [inner_radius,height];
% initial_state{4} = [inner_radius,height];
% initial_state{5} = [outter_radius,height];
% initial_state{6} = [outter_radius,height];
% initial_state{7} = [outter_radius,height];
% initial_state{8} = [outter_radius,height];