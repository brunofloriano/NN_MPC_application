function VideoName = plot_movie(savefolderlocal,MP4,GIF,t,individual_coord,Ggraph,GRAPH_AREA,Nr,encompassed_area_time)

if GRAPH_AREA == 1
    video_file_name = '/HurricaneVideoArea2.mp4';
else
    video_file_name = '/HurricaneVideoGraph.mp4';
end

%file_name = [save_mainname time_now]; %'dataBalloon 2022-12-13-11-50';
%load(['results/' file_name '/' file_name '.mat']);
mkdir([savefolderlocal '/frames_graph']);
if MP4 == 1
    VideoName = VideoWriter([savefolderlocal video_file_name],'MPEG-4');
    open(VideoName);
end
% for agent_counter = 1:N
%     plot_window{agent_counter} = plot(NaN,NaN);
% end
figure
%hold on
grid
%axis(200e3*[-1 1 -1 1]);

%list = 'rbkgbmyrb';

n_steps = 1;%length(t)/1000;%200;
counter = 0;
for time_counter = 1:n_steps:length(t)
    counter = counter + 1;

    if GRAPH_AREA == 1
        mesh(encompassed_area_time{time_counter});
        %contour(encompassed_area_time{time_counter},Nr+1)
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
    end
    title([ 't = ' num2str( t(time_counter) ) ' seconds'])
    print([savefolderlocal '/frames_graph/Frame ' num2str(counter)], '-dpng', '-r150');

    if MP4 == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
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
if MP4 == 1
    close(VideoName);
end