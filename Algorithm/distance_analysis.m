clear all; close all; clc;

%file_name = 'dataBalloon 2022-10-20-14-47';
%file_name = 'dataBalloon 2022-11-1-10-15';
file_name = 'dataBalloon 2022-11-4-16-37';
load(['results/' file_name '/' file_name '.mat']);

number_pairs = nchoosek(N,2);
%distance = zeros(number_pairs,1);
list = 'rbkgbmyr';

for time_counter = 1:length(t)
    pair_counter = 0;
    for agent_count_i = 5%1:N
        for agent_count_j = 1:N
            if agent_count_i ~= agent_count_j%agent_count_i < agent_count_j
                pair_counter = pair_counter + 1;
                distance{pair_counter}(time_counter) = norm( individual_coord{agent_count_i}(:,time_counter) - individual_coord{agent_count_j}(:,time_counter) );
            end
        end
    end
end
number_pairs = pair_counter;
figure; grid; hold on;
for pair_counter = 1:number_pairs
    plot(t,distance{pair_counter}/1e3,'LineWidth',2)
end

xlabel('Time (s)')
ylabel('Distance (km)')
xlim([0 t(end)])

% saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_distances.eps']);
% saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_distances.fig']);
saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_distances.jpg']);


% if time_counter > 1
%     velocity{pair_counter}(time_counter) = ( distance{pair_counter}(time_counter) - distance{pair_counter}(time_counter-1) )/delta_t;
% else
%     velocity{pair_counter}(time_counter) = 0;
% end

shift_size = 5000;
figure; grid; hold on;

for pair_counter = 1:number_pairs
    vector = distance{pair_counter};

    shifted_vector_begin = [zeros(1,shift_size) vector];
    shifted_vector_end = [vector zeros(1,shift_size)];

    difference_vector = shifted_vector_end - shifted_vector_begin;
    velocity{pair_counter} = [zeros(1,shift_size) difference_vector(shift_size+1:end-shift_size)]/(delta_t*shift_size);

    plot(t,velocity{pair_counter}*3.6);
end

xlabel('Time (s)')
ylabel('Velocity (km/h)')
xlim([0 t(end)])

% saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_velocities.eps']);
% saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_velocities.fig']);
saveas(gcf,[savefolderlocal '/' save_mainname  time_now '_velocities2.jpg']);
