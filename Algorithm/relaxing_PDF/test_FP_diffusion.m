% Clear workspace
clear; clc; close all;

% Initiate video parameters
VIDEO = 1;
if VIDEO == 1
    clocktime = clock;
    time_now = [num2str(clocktime(1)) '-' num2str(clocktime(2)) '-' num2str(clocktime(3)) '-' num2str(clocktime(4)) '-' num2str(clocktime(5))];
    video_file_name = ['DiffusionTest' time_now '.mp4'];
    VideoName = VideoWriter([video_file_name],'MPEG-4');
    open(VideoName);
end

% Add functions
addpath functions

% Initialize time parameters
delta_t = 0.1;
tmax = 1000;
t = 0:delta_t:tmax;

% Initialize space parameters
delta_x = 0.1;
xmax = 15;
xmaxvisual = 10;
x = -xmax:delta_x:xmax;
xvisual = -xmaxvisual:delta_x:xmaxvisual;

%% Initialize PDF as a uniform distribution

p_initial = 0*x;
counter = 0;
for xcounter = 1:length(x)
    if x(xcounter) > - 5 && x(xcounter) < 5
        p_initial(xcounter) = 1/10;
        counter = counter + 1;
    else
        p_initial(xcounter) = 0;
    end
end
%p_initial = p_initial/counter;
%p_initial = normpdf(x,0,5);

p_desired = normpdf(x,0,1);


% No drift
D = 1e-2;
u = -D*x.^3;
%u(abs(x)>20) = 0;

%% Time loop
for time_counter = 1:length(t)
    if time_counter > 1
        %D = 0.5;
        %D = 0.5*(10 - std(p{time_counter-1}))/(10 - std(p{1}));
        %D = 2*norm(p{time_counter-1} - p_desired);
        %D = 2*(p{time_counter-1} - p_desired);

        p{time_counter} = FP_diffusion(p{time_counter-1},D,u,delta_t,delta_x);
    else
        p{time_counter} = p_initial;
    end
%     p{time_counter}(1) = 0;
%     p{time_counter}(end) = 0;
    p{time_counter}(x>xmaxvisual) = 0;
    p{time_counter}(x<-xmaxvisual) = 0;
    
    % p{time_counter} = p{time_counter}/max(p{time_counter});
    plot(x,p{time_counter})
    hold on
    plot(x,p_desired);
    xlim([-xmaxvisual xmaxvisual])
    ylim([0 0.5])
    hold off
    if VIDEO == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
    end
end
if VIDEO == 1
    close(VideoName);
end