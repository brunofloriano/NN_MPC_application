clear all; close all; clc;

load('results/dataBalloon 2022-10-20-14-37/dataBalloon 2022-10-20-14-37.mat');


plot_window = plot(NaN,NaN);
hold on
axis(90e3*[-1 1 -1 1]);
list = 'rbkgwmyc';

for time_counter = 1:200:length(t)
    %pause(0.01)
    for agent_counter = 1:N
        x = individual_coord{agent_counter}(1,1:time_counter);
        y = individual_coord{agent_counter}(2,1:time_counter);
        
        %plot(x,y,list(agent_counter))
        set(plot_window, 'XData', x, 'YData', y, 'LineWidth',3, 'Color',list(agent_counter));
        drawnow
    end
end