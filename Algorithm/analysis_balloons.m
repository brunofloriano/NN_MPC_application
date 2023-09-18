clear; close all; clc;
addpath 'functions'
addpath 'results'

%file_name = 'dataBalloon 2023-4-24-22-45';
%file_name = 'dataBalloon 2023-4-27-3-15'; % 24h, 500 steps, D=1
%file_name = 'dataBalloon 2023-5-6-4-6'; % 24h + 24h, 500 + 500 steps, D=1
%file_name = 'dataBalloon 2023-5-8-17-46'; % 24h + 24h, 500 + 500 steps, D=1, zbar corrected
%file_name = 'dataBalloon 2023-5-9-11-9'; % 48h, 1000 steps, D=0.5
%file_name = 'dataBalloon 2023-5-11-6-33'; % 48h, 500 steps, D=1
%file_name = 'dataBalloon 2023-5-11-18-28'; % 48h, 500 steps, D=1, communication decay, tf=12h
%file_name = 'dataBalloon 2023-5-12-5-17'; % 48h, 500 steps, D=1, communication decay, tf=6h
%file_name = 'dataBalloon 2023-5-12-10-0'; % 48h, 50 steps, D=1, communication decay corrected, tf=6h

file_name = 'dataBalloon 2023-5-12-18-51'; cycles = 48/6; interval = 75; ROUTING = 0; CUMMULATIVE = 1; % 48h, 500 steps, D=1, communication decay corrected, tf=6h
%file_name = 'dataBalloon 2023-5-13-6-17'; cycles = 48/12; interval = 125; ROUTING = 0; CUMMULATIVE = 1; % 48h, 500 steps, D=1, communication decay corrected, tf=12h
%file_name = 'dataBalloon 2023-5-30-23-16'; cycles = 48/6; interval = 75; ROUTING = 1; CUMMULATIVE = 0; % N=3+1, 48h, 500 steps, D=1, communication decay corrected, tf=6h

file = ['results/' file_name '/' file_name '.mat'];

load(file);

AREACOMM = 0;
DISTANCE = 0;
RADIUS = 0;
%CUMMULATIVE = 0;
TABLE = 0;
%ROUTING = 0;
POSITION = 0;
POSITION2 = 0;
DISTANCE2 = 0;

%% Area covered with communication
if AREACOMM == 1
    figure

    for time_counter = 1:length(t)
        parameter(time_counter) = 1 - sum(sum(encompassed_area_time{time_counter}))/sum(sum(Area.initial_encompassed_area));

    end
    parameter = parameter + abs(min(parameter)) + 0.05;
    %parameter = (parameter)/max(parameter);

    yyaxis left
    plot(t/(3600),parameter)
    grid
    xlabel('t')
    ylabel('A_c','FontName', 'Arial')

    for time_counter = 1:length(t)
        parameter_cc(time_counter) = -beta_cc * ( L_cumulative_time{time_counter}(1,1) + gamma_cc*( sum(diag(L_cumulative_time{time_counter})) -  L_cumulative_time{time_counter}(1,1)) ); %L_cumulative_time{time_counter}(1,1);
    end

    yyaxis right
    plot(t/(3600),parameter_cc,'--')
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$\bar{l}_{cc}$','FontName', 'Arial','Interpreter','Latex')
    %ylim([-0.5 8.5])

    name = 'AreaComm';
    save_figure(gcf,name,save_folder,save_mainname,time_now);
    
    counter = 1;
    for time_counter = [1 100 200 300 400 501]
        Vac(counter) = 1 - parameter(time_counter);
        Vcc(counter) = parameter_cc(time_counter);
        vector_norm = norm([Vac(counter) Vcc(counter)]);
        importance_size(counter) = (Vcc(counter))/vector_norm
        counter = counter + 1;
    end
end

%% Cummulative area coverage
if CUMMULATIVE == 1
    %cycles = Time.tmax/(6*3600);
    %interval = round(length(t)/cycles);%75;%125;%167;%round(length(t)/cycles);
    total_sum = sum(sum(Area.initial_encompassed_area));
    counter = 0;
    for time_counter = 1:length(t)
        if counter >= interval
            %counter = 0;
            condition2 = cummulative_area < cummulative_area_time{time_counter-interval+1};
            cummulative_area(condition2) = cummulative_area_time{time_counter-interval+1}(condition2);
        end

        if counter == 0
            cummulative_area = encompassed_area_time{time_counter};
        else
            condition = cummulative_area > encompassed_area_time{time_counter};
            cummulative_area(condition) = encompassed_area_time{time_counter}(condition);
        end
        cummulative_area_time{time_counter} = cummulative_area;
        counter = counter + 1;

        An_bar(time_counter) = sum(sum(cummulative_area_time{time_counter}));
        Ac_bar(time_counter) =  1 - An_bar(time_counter)/total_sum;
    end

    figure
    yyaxis left
    plot(t/(3600),Ac_bar)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$\bar{A}_c$','FontName', 'Arial','Interpreter','latex')
    %ylim([0.2 0.55])

    for time_counter = 1:length(t)
        parameter_cc(time_counter) = -beta_cc * ( L_cumulative_time{time_counter}(1,1) + gamma_cc*( sum(diag(L_cumulative_time{time_counter})) -  L_cumulative_time{time_counter}(1,1)) ); %L_cumulative_time{time_counter}(1,1);
    end

    yyaxis right
    plot(t/(3600),parameter_cc,'--')
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$E_{cc}$','FontName', 'Arial','Interpreter','Latex')
    %ylim([-0.5 8.5])
    %fontsize(16,"points")

    name = 'CummulativeAreaComm';
    save_figure(gcf,name,save_folder,save_mainname,time_now);

    %     interval = 10; summ = 0; total_sum = sum(sum(Area.initial_encompassed_area));
    %     for time_counter = 1:length(t)
    %         An(time_counter) = sum(sum(encompassed_area_time{time_counter}));
    %         if time_counter <= interval
    %             summ = summ + An(time_counter)*delta_t;
    %         else
    %             summ = summ + An(time_counter)*delta_t - An(time_counter-interval)*delta_t;
    %         end
    %         An_bar(time_counter) = summ/(interval*delta_t);
    %         Ac_bar(time_counter) =  1 - An_bar(time_counter)/total_sum;
    %     end
    %
    %     figure
    %     hold on
    %     plot(t/(3600),Ac_bar/max(Ac_bar))
    %     grid
    %     xlabel('t(h)','FontName', 'Arial')
    %     ylabel('$\bar{A}_c$','FontName', 'Arial','Interpreter','latex')
end


%% Area covered over time
% figure
%
% for time_counter = 1:length(t)
%     parameter(time_counter) = 1 - sum(sum(encompassed_area_time{time_counter}))/sum(sum(Area.initial_encompassed_area));
%
% end
% parameter = parameter + abs(min(parameter)) + 0.05;
% %parameter = (parameter)/max(parameter);
%
% plot(t/(3600),parameter)
% grid
% xlabel('t')
% ylabel('Eac')

%% Communication
% figure
%
% for time_counter = 1:length(t)
%     parameter(time_counter) = L_cumulative_time{time_counter}(1,1);
% end
%
% plot(t/(3600),parameter)
% grid
% xlabel('t')
% ylabel('lcc')
% ylim([-0.5 7.5])

%% Position
if POSITION == 1
id7 = 1; id8 = 1; id1 = 1;
for i = 1:N
    for time_counter = 1:length(t)
        x = individual_coord{i}(1,time_counter);
        y = individual_coord{i}(2,time_counter);

        r(time_counter,i) = norm([x y]);
        theta(time_counter,i) = atan2(y,x);

        x1 = individual_coord{1}(1,time_counter);
        y1 = individual_coord{1}(2,time_counter);
        di1(i,time_counter) = norm([x y] - [x1 y1]);

        x7 = individual_coord{7}(1,time_counter);
        y7 = individual_coord{7}(2,time_counter);
        di7(i,time_counter) = norm([x y] - [x7 y7]);


        x8 = individual_coord{8}(1,time_counter);
        y8 = individual_coord{8}(2,time_counter);
        di8(i,time_counter) = norm([x y] - [x8 y8]);


    end
    if min(di7(i,:)) < 10e3
        include7(id7) = i;
        id7 = id7 + 1;
    end

    if min(di8(i,:)) < 10e3
        include8(id8) = i;
        id8 = id8 + 1;
    end

    if min(di1(i,:)) < 10e3
        include1(id1) = i;
        id1 = id1 + 1;
    end
end
end
% Distance plot

if DISTANCE == 1
    figure
    hold on
    for i = 5
        plot(t/(3600),di7(i,:)/1e3);
        plot(t/(3600),di1(i,:)/1e3);
    end
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('Distance (km)','FontName', 'Arial')
    legend('d57','d15','FontName', 'Arial')

    figure
    hold on
    for i = 2
        plot(t/(3600),di8(i,:)/1e3);
        plot(t/(3600),di1(i,:)/1e3);
    end
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('Distance (km)','FontName', 'Arial')
    legend('d28','d25','FontName', 'Arial')

    name = 'Distance';
    save_figure(gcf,name,save_folder,save_mainname,time_now);

end

% % Single plot
%
% figure
% hold on
% for agent_counter = 2:N
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
%
% % Double plot
% figure
% subplot(2,1,1)
% hold on
% for agent_counter = 2:5
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
% grid
% xlabel('t')
% ylabel('r')
% subplot(2,1,2)
% hold on
% for agent_counter = 6:9
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
% grid
% xlabel('t')
% ylabel('r')
%
% % Triple plot
% figure
% subplot(3,1,1)
% hold on
% for agent_counter = 2:4
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
% plot(t/(3600),t*0 + 10)
% grid
% xlabel('t')
% ylabel('r')
% subplot(3,1,2)
% hold on
% for agent_counter = 5:7
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
% plot(t/(3600),t*0 + 10)
% grid
% xlabel('t')
% ylabel('r')
% legend('1','2','3')
% subplot(3,1,3)
% hold on
% for agent_counter = 8:9
%     plot(t/(3600),r(:,agent_counter)/1e3);
% end
% plot(t/(3600),t*0 + 10)
% grid
% xlabel('t')
% ylabel('r')

% % Differenticate agents plot
% legends{1} = '-';
% legends{2} = '--';
% legends{3} = ':';
% legends{4} = '-.';
% legends{5} = '-';
% figure
% subplot(2,1,1)
% hold on
% for agent_counter = 2:6
%     plot(t/(3600),r(:,agent_counter)/1e3,legends{agent_counter-1});
% end
% plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
% grid
% xlabel('t')
% ylabel('r')
% legend('a2','a3','a5','a5','a6')
% subplot(2,1,2)
% hold on
% for agent_counter = 7:9
%     plot(t/(3600),r(:,agent_counter)/1e3,legends{agent_counter-6});
% end
% plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
% grid
% xlabel('t')
% ylabel('r')
% legend('a7','a8','a9')

% Differenticate agents three plot
if RADIUS == 1
    legends{1} = '-';
    legends{2} = '--';
    legends{3} = '-.';
    legends{4} = ':';
    legends{5} = '-';
    figure

    subplot(3,1,1)
    hold on
    for i = 2:4
        plot(t/(3600),r(:,i)/1e3,legends{i-1});
    end
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('Radius (km)','FontName', 'Arial')
    legend('a2','a3','a4','FontName', 'Arial')
    yticks([10 50 100 150])

    subplot(3,1,2)
    hold on
    for i = 5:6
        plot(t/(3600),r(:,i)/1e3,legends{i-4});
    end
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('Radius (km)','FontName', 'Arial')
    legend('a5','a6','FontName', 'Arial')
    yticks([10 50 100 150 200])


    subplot(3,1,3)
    hold on
    for i = 7:9
        plot(t/(3600),r(:,i)/1e3,legends{i-6});
    end
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('Radius (km)','FontName', 'Arial')
    legend('a7','a8','a9','FontName', 'Arial')
    yticks([10 50 100 150 200 250])

    name = 'Radius';
    save_figure(gcf,name,save_folder,save_mainname,time_now);
end

%Angle
% figure
% hold on
% for agent_counter = 1:N
%     plot(t,theta(:,agent_counter));
% end
% grid
% xlabel('t')
% ylabel('theta')

%% Communication table

if TABLE == 1
    ncolumns = 4;
    %tablenum = zeros(N,ncolumns);
    %tablechar = num2str(zeros(N,ncolumns*N));
    %tabelchar = [];
    step = round(length(t)/ncolumns);
    counter = 0;
    %for time_counter = 1:length(t)-2
        %counter = floor(time_counter/step) + 1;
    for time_counter = [1 125 250 375 500]
        counter = counter + 1;

        for i = 1:N
            word = []; connect = 0;
            for j = 2:N
                item = L_cumulative_time{time_counter}(i,j);
                if item < 0
                    word = [word num2str(j) ' '];
                    connect = 1;
                end
            end
            if connect == 0
                word = 'Empty';
            end
            table{i,counter} = word;
        end
    end
    table
end


% Differenticate agents three plot
if ROUTING == 1
    legends{1} = '-';
    legends{2} = '--';
    legends{3} = '-.';
    legends{4} = ':';
    legends{5} = '-';
    figure

    parameter_l = zeros(length(t),N,N);
    for time_counter = 1:length(t)
        parameter_l(time_counter,:,:) = L_cumulative_time{time_counter};
        %parameter_l(time_counter,:,:) = Laplacian_time{time_counter};
    end

    subplot(2,1,1)
    hold on
    for i = 4:4
        plot(t/(3600),parameter_l(:,1,i),legends{i-1});
    end
    %plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$\bar{l}_{ij}$','FontName', 'Arial','Interpreter','Latex')
    legend('$\bar{l}_{14}$','FontName', 'Arial','Interpreter','Latex')
    %yticks([10 50 100 150])

    subplot(2,1,2)
    hold on
    %plot(t/(3600),parameter_l(:,2,3),legends{1});
    plot(t/(3600),parameter_l(:,2,4),legends{2});
    plot(t/(3600),parameter_l(:,3,4),legends{3});

    %plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$\bar{l}_{ij}$','FontName', 'Arial','Interpreter','Latex')
    legend('$\bar{l}_{24}$','$\bar{l}_{34}$','FontName', 'Arial','Interpreter','Latex')
    %yticks([10 50 100 150 200])


%     subplot(3,1,3)
%     hold on
%     for agent_counter = 7:9
%         plot(t/(3600),r(:,agent_counter)/1e3,legends{agent_counter-6});
%     end
%     plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
%     grid
%     xlabel('t(h)','FontName', 'Arial')
%     ylabel('Radius (km)','FontName', 'Arial')
%     legend('a7','a8','a9','FontName', 'Arial')
%     yticks([10 50 100 150 200 250])

    name = 'Routing';
    save_figure(gcf,name,save_folder,save_mainname,time_now);
end

%% Position
if POSITION2 == 1
    dij = zeros(length(t),N,N);
%for agent_counter = 1:N
    for time_counter = 1:length(t)
        for i = 1:N
            x(i) = individual_coord{i}(1,time_counter);
            y(i) = individual_coord{i}(2,time_counter);
        end

        for i = 1:N
            for j = 1:N
                 dij(time_counter,i,j) = norm([x(i) y(i)] - [x(j) y(j)]);
            end
        end
    end
%end 
end
% Distance plot

if DISTANCE2 == 1
    figure
    hold on
    plot(t/(3600),dij(:,1,2)/1e3,legends{1});
    plot(t/(3600),dij(:,1,3)/1e3,legends{2});
    plot(t/(3600),dij(:,1,4)/1e3,legends{3});
    plot(t/(3600),dij(:,2,3)/1e3,legends{4});
    plot(t/(3600),dij(:,2,4)/1e3,legends{5});
    plot(t/(3600),dij(:,3,4)/1e3,legends{1});
    plot(t/(3600),t*0 + 10,'k','LineWidth',1.5)
    grid
    xlabel('t(h)','FontName', 'Arial')
    ylabel('$d_{ij}$','FontName', 'Arial','Interpreter','Latex')
    legend('$\bar{d}_{12}$','$\bar{d}_{13}$','$\bar{d}_{14}$','$\bar{d}_{23}$','$\bar{d}_{24}$','$\bar{d}_{34}$','FontName', 'Arial','Interpreter','Latex')
    %yticks([10 50 100 150])

    name = 'Distance2';
    save_figure(gcf,name,save_folder,save_mainname,time_now);

end
