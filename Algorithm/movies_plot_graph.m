clear all; close all; clc;

%file_name = 'dataBalloon 2022-10-20-14-47';
%file_name = 'dataBalloon 2022-11-1-10-15';
%file_name = 'dataBalloon 2022-11-15-16-4';
%file_name = 'dataBalloon 2022-11-28-14-51';
file_name = 'dataBalloon 2022-11-28-15-37';
load(['results/' file_name '/' file_name '.mat']);
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
    title([ 't = ' num2str( t(time_counter) ) ' seconds'])
    print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');

    disp(['Simulation at ' num2str(time_counter*100/length(t)) ' %' newline])

end

%% Now time to generate the animated gif
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