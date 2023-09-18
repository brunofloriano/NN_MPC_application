function VideoName = plot_movie(savefolderlocal,VIDEO,GIF,t,individual_coord,Ggraph,GRAPH_AREA,Nr,encompassed_area_time)

if GRAPH_AREA == 1
    %video_file_name = '/HurricaneVideoArea2.avi';
    video_file_name = '/HurricaneVideoArea1-2.mp4';
    video_file_name2 = '/HurricaneVideoArea2-2.mp4';
else
    video_file_name = '/HurricaneVideoGraph2.avi';
end

%file_name = [save_mainname time_now]; %'dataBalloon 2022-12-13-11-50';
%load(['results/' file_name '/' file_name '.mat']);
mkdir([savefolderlocal '/frames_graph']);
mkdir([savefolderlocal '/frames_graph/FramesOut']);
if VIDEO == 1
    %VideoName = VideoWriter([savefolderlocal video_file_name],'Motion JPEG AVI');
    VideoName = VideoWriter([savefolderlocal video_file_name],'MPEG-4');
    VideoName2 = VideoWriter([savefolderlocal video_file_name2],'MPEG-4');
    open(VideoName);
    open(VideoName2);
end
% for agent_counter = 1:N
%     plot_window{agent_counter} = plot(NaN,NaN);
% end
%figure
%hold on
%grid
%axis(200e3*[-1 1 -1 1]);

%list = 'rbkgbmyrb';

n_steps = 1;%length(t)/1000;%200;
counter = 0;
f1 = figure;
for time_counter = 1:n_steps:length(t)
    counter = counter + 1;

    if GRAPH_AREA == 1
        %% Imagesc plot

        %mesh(encompassed_area_time{time_counter});
        %contour(encompassed_area_time{time_counter},Nr+1)
        %contourf(encompassed_area_time{time_counter})
        imagesc([-200e3 200e3]/1e3, [-200e3 200e3]/1e3, encompassed_area_time{time_counter})
        %imagesc(encompassed_area_time{time_counter})
        %heatmap(encompassed_area_time{time_counter}); grid off;

        hold on

        for agent_counter = 1:Nr+1
            x = individual_coord{agent_counter}(1,time_counter);
            y = individual_coord{agent_counter}(2,time_counter);
            graph_position(:,agent_counter) = [x;y];
            %plot(x,y, 'LineWidth',3)

            %set(plot_window, 'XData', x, 'YData', y, 'LineWidth',3, 'Color',list(agent_counter));
            %drawnow
        end

        % NO BALLOON NUMBER
        %plot(Ggraph{time_counter},'NodeLabel',{},'NodeColor','k',XData=graph_position(1,:)/1e3,YData=real(graph_position(2,:))/1e3)
        % BALLOON NUMBER
        plot(Ggraph{time_counter},'NodeColor','k',XData=graph_position(1,:)/1e3,YData=real(graph_position(2,:))/1e3)
        
        
        %fontsize(gcf,64,"points")
        %fontsize(gcf,scale=1.5)
        %axis(200e3*[-1 1 -1 1]);
        hold off

        %title([ 't = ' num2str( t(time_counter) ) ' seconds'])

        title([ 't = ' num2str( t(time_counter)/3600 ) ' hours'])
        xlabel('X (km)')
        ylabel('Y (km)')

        % UNCOMMENT HERE TO SET XTICK AND YTICK EMPTY
        % set(gca,'xtick',[])
        % set(gca,'ytick',[])

        %print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');
        saveas(gcf,[savefolderlocal '/frames_graph/FramesOut/Frame_' num2str(counter)],'epsc')

        
        if VIDEO == 1
            VideoFrame = getframe(f1);
            writeVideo(VideoName,VideoFrame);
        end
        if time_counter == 245
            %saveas(gcf,[savefolderlocal '/frames_graph/communication2d.jpg'])
            saveas(gcf,[savefolderlocal '/frames_graph/communication2d'],'epsc')
        end

    else
        %pause(0.01)
        for agent_counter = 1:Nr+1
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
        %print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');

        if VIDEO == 1
            VideoFrame = getframe(gcf);
            writeVideo(VideoName,VideoFrame);
        end
    end

    disp(['Graph plotting progress: ' num2str(time_counter*100/length(t)) ' %' newline])

end

n_steps = 1;%length(t)/1000;%200;
counter = 0;
f2 = figure;
for time_counter = 1:n_steps:length(t)
    counter = counter + 1;

    if GRAPH_AREA == 1
        %% Mesh plot

        mesh(encompassed_area_time{time_counter});

        %title([ 't = ' num2str( t(time_counter) ) ' seconds'])
        title([ 't = ' num2str( t(time_counter)/3600 ) ' hours'])
        %print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');

        if VIDEO == 1
            VideoFrame = getframe(f2);
            writeVideo(VideoName2,VideoFrame);
        end
    end

    disp(['Graph plotting progress: ' num2str(time_counter*100/length(t)) ' %' newline])

end

if GIF == 1
    % Now time to generate the animated gif
    GifName = [savefolderlocal '/HurricaneAnimation.gif'];

    delay = 0.1;    % Delay between frames (s)
    for ii = 1:1000
        [A, ~] = imread([savefolderlocal '/frames_graph/Frame ' num2str(ii) '.png']);
        [X, map] = rgb2ind(A, 256);
        if ii == 1
            imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
        else
            imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
        end
    end
end
if VIDEO == 1
    close(VideoName);
    if GRAPH_AREA == 1
        close(VideoName2);
    end
end