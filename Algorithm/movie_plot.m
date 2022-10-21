clear all; close all; clc;

load('results/dataBalloon 2022-10-20-14-37/dataBalloon 2022-10-20-14-37.mat');

% for agent_counter = 1:N
%     plot_window{agent_counter} = plot(NaN,NaN);
% end
figure
hold on
grid
axis(90e3*[-1 1 -1 1]);
list = 'rbkgbmyr';

n_steps = length(t)/200;
counter = 0;
for time_counter = 1:n_steps:length(t)
    counter = counter + 1;
    %pause(0.01)
    for agent_counter = 1:N
        x = individual_coord{agent_counter}(1,1:time_counter);
        y = individual_coord{agent_counter}(2,1:time_counter);
        
        plot(x,y, 'LineWidth',3, 'Color',list(agent_counter))
        title([ 't = ' num2str( t(time_counter) ) ' seconds'])
        %set(plot_window, 'XData', x, 'YData', y, 'LineWidth',3, 'Color',list(agent_counter));
        %drawnow
    end
    print(['frames/Frame ' num2str(counter)], '-dpng', '-r150');
end

%% Now time to generate the animated gif
GifName = 'HurricaneAnimation.gif';
delay = 0.1;    % Delay between frames (s)
for ii = 1:n_steps
    [A, ~] = imread(['frames/Frame ' num2str(ii) '.png']);
    [X, map] = rgb2ind(A, 256);
    if ii == 1
        imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
    else
        imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
    end
end