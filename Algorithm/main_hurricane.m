%% This code implements a group of balloons around a hurricane
%% by Bruno R. O. Floriano

clear; close all; clc;
addpath 'functions'

%% Include the system functions and parameters
addpath 'Balloon'
start
N = Nr;
save_folder = 'results/';

%% Choose simulation parameters
SAVE_DATA = 1;
GRAPH_MOVIE = 1;

%% Generate initial conditions (inner and outer radius)
% N = length(initial_conditions);
% 
% initial_angular_position = initial_angular_positions();
% initial_condition = initial_conditions();
% 
% individual_angular_position = initial_angular_position;
% individual_state = initial_condition;
% 
% for agent_counter = 1:N
%     [x,y] = pol2cart(individual_angular_position{agent_counter},individual_state{agent_counter}(1));
%     individual_coord{agent_counter}(:,1) = [x;y];
% end

%% Generate initial conditions (random agents)

agents_position = generate_agents(N,max_radius);
position = [[0;0] agents_position];
[agents_angular_position,agents_radius] = cart2pol(position(1,:),position(2,:));

for agent_counter = 1:N+1
    individual_coord{agent_counter}(:,1) = position(:,agent_counter);
    individual_angular_position{agent_counter} = agents_angular_position(agent_counter);
    individual_state{agent_counter} = [agents_radius(agent_counter);height];
end

%% Simulation
t = 0:delta_t:tmax-delta_t;

Ggraph{1} = routing_protocol(position,comm_range,max_root_connections);
% Time loop
for time_counter = 2:length(t)
    for agent_counter = 1:N+1
        control_input = 0;
        if individual_state{agent_counter}(1) > 20e3
            control_input = -level_height;
        end
        individual_state{agent_counter} = dynamics(individual_state{agent_counter},control_input,delta_t,level_height);
        
        theta = individual_angular_position{agent_counter};
        individual_angular_vector{agent_counter}(time_counter) = theta;

        radius = individual_state{agent_counter}(1);
        
        if radius < 0
            %theta_position = - theta_position;
            radius = - radius;
            %individual_state{agent_counter+1}(1) = - individual_state{agent_counter+1}(1);
            %disp('Ta dando treta aqui')
        end
        radius_vector{agent_counter}(time_counter) = radius;

        [x,y,theta,radius] = kinematics(theta,radius,delta_t);
        x = real(x);
        y = real(y);
        position(:,agent_counter) = [x;y];
        [theta_position,radius] = cart2pol(x,y);

        individual_coord{agent_counter}(:,time_counter) = [x;y];
        individual_angular_position{agent_counter} = theta;
    end
    Ggraph{time_counter} = routing_protocol(position,comm_range,max_root_connections);

    disp(['Simulation at ' num2str(time_counter*100/length(t)) ' %' newline])
end

% for agent_counter = 1:N
%     plot(individual_coord{agent_counter}(1,:),individual_coord{agent_counter}(2,:))
%     hold on
% end


disp(['Saving data' newline])
%% Save data
if SAVE_DATA == 1 || GRAPH_MOVIE == 1
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

disp(['Generating graph plot'])
%% Generate movie plot
if GRAPH_MOVIE == 1
    %file_name = [save_mainname time_now]; %'dataBalloon 2022-12-13-11-50';
    %load(['results/' file_name '/' file_name '.mat']);
    mkdir([savefolderlocal '/frames_graph']);

    % for agent_counter = 1:N
    %     plot_window{agent_counter} = plot(NaN,NaN);
    % end
    figure
    %hold on
    grid
    axis(200e3*[-1 1 -1 1]);

    %list = 'rbkgbmyrb';

    n_steps = length(t)/200;
    counter = 0;
    for time_counter = 1:n_steps:length(t)
        counter = counter + 1;
        %pause(0.01)
        for agent_counter = 1:N+1
            x = individual_coord{agent_counter}(1,time_counter);
            y = individual_coord{agent_counter}(2,time_counter);
            graph_position(:,agent_counter) = [x;y];
            %plot(x,y, 'LineWidth',3)

            %set(plot_window, 'XData', x, 'YData', y, 'LineWidth',3, 'Color',list(agent_counter));
            %drawnow
        end
        plot(Ggraph{time_counter},XData=graph_position(1,:),YData=real(graph_position(2,:)))
        axis(200e3*[-1 1 -1 1]);
        title([ 't = ' num2str( t(time_counter) ) ' seconds'])
        print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');

        disp(['Graph plotting progress: ' num2str(time_counter*100/length(t)) ' %' newline])

    end

    % Now time to generate the animated gif
    GifName = [savefolderlocal '/HurricaneAnimation.gif'];
    delay = 0.1;    % Delay between frames (s)
    for ii = 1:n_steps
        [A, ~] = imread([savefolderlocal '/frames_graph/Frame ' num2str(ii) '.png']);
        [X, map] = rgb2ind(A, 256);
        if ii == 1
            imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
        else
            imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
        end
    end
end