% Clear workspace
clear; clc; close all;

% Initiate video parameters
VIDEO = 1;
SAVE_VIDEO = 0;

% Video parameters
if VIDEO == 1
    clocktime = clock;
    if SAVE_VIDEO == 1
        time_now = [num2str(clocktime(1)) '-' num2str(clocktime(2)) '-' num2str(clocktime(3)) '-' num2str(clocktime(4)) '-' num2str(clocktime(5))];
        video_file_name = ['DiffusionTest2D ' time_now '.mp4'];
    else
        video_file_name = ['DiffusionTest2D.mp4'];
    end
    VideoName = VideoWriter([video_file_name],'MPEG-4');
    open(VideoName);
end

% Add functions
addpath functions

% Initialize time parameters
delta_t = 0.001;
tmax = 10;
t = 0:delta_t:tmax;

% Initialize space parameters
delta_x = 0.1;
xmax = 15;
xmaxvisual = xmax;
x = -xmax:delta_x:xmax;
y = x;
xvisual = -xmaxvisual:delta_x:xmaxvisual;
yvisual = xvisual;
epsilon = 1e-6;

[x1,x2] = meshgrid(x,y);
xy = [x1(:) x2(:)];
n = length(x);

noboundary = ones(n);
square_boundary = ones(n);
square_boundary(:,1) = zeros(n,1);
square_boundary(:,end) = zeros(n,1);
square_boundary(1,:) = zeros(1,n);
square_boundary(end,:) = zeros(1,n);

Dn = 1e-0;
D = Dn*[1 0;0 1];
%D = Dn*[1 1;1 1];

%% Initialize PDF as a uniform distribution

p_initial = zeros(n);
inside = p_initial;
round_boundary = inside;
un = p_initial;
uniform_lim = 5;
counter = 0;
for xcounter = 1:n
    for ycounter = 1:n
        xvector = [x(xcounter); y(ycounter)];
        r = norm(xvector);
        theta = atan2(y(ycounter),x(xcounter));

        p_initial(xcounter,ycounter) = (1 - 1./(1+exp(-(r-0/2))))*4/100;
        
                if abs(x(xcounter)) < uniform_lim && abs(y(ycounter)) < uniform_lim
                    p_initial(xcounter,ycounter) = 1/(uniform_lim*2)^2;
                    counter = counter + 1;
                else
                    p_initial(xcounter,ycounter) = 0;
                end

        %un(xcounter,ycounter,1:2) = -D*xvector;
        %un(xcounter,ycounter,1:2) = -D*xvector.^3;
        un(xcounter,ycounter,1:2) = -D*r^3*[cos(theta); sin(theta)];
        %un(xcounter,ycounter,1:2) = -D*r^2*xvector;

        if abs(x(xcounter)) < xmaxvisual && abs(y(ycounter)) < xmaxvisual
            inside(xcounter,ycounter) = 1;
        end
        if r < xmax - delta_x*0
            round_boundary(xcounter,ycounter) = 1;
        end

    end
end

p_desired = mvnpdf(xy,[0 0],[1 0;0 1]);
p_desired = reshape(p_desired,length(x),length(x));

p_initial = mvnpdf(xy,[0 0],[.5 0;0 .5]);
p_initial = reshape(p_initial,length(x),length(x));

% No drift

u = un;

%boundary = noboundary;
boundary = square_boundary;
%boundary = round_boundary;
%boundary = inside;

%% Time loop
%for time_counter = 1:length(t)
time_counter = 0;
current_epsilon = 0;
old_epsilon = 1;
%while true %current_epsilon < old_epsilon
%time_counter = time_counter + 1;
for time_counter = 1:length(t)
    %boundary(:,time_counter) = zeros(n,1);
    %boundary(:,n-time_counter+1) = zeros(n,1);
    %boundary(time_counter,:) = zeros(1,n);
    %boundary(n-time_counter+1,:) = zeros(1,n);
    if time_counter > 1

        %p{time_counter} = FP_diffusion2D(p{time_counter-1},D,u,delta_t,delta_x,p_initial,boundary);
        p{time_counter} = FP_diffusion2D2(p{time_counter-1},D,u,delta_t,delta_x,boundary);
        pd{time_counter} = p{time_counter} - p{time_counter-1};

        current_epsilon = max(max(pd{time_counter}));
        old_epsilon = max(max(pd{time_counter - 1}));
    else
        p{time_counter} = p_initial;
        pd{time_counter} = p{time_counter};
    end


    p{time_counter} = p{time_counter}.*boundary;

    mesh(x,y,p{time_counter})
    %imagesc(xvisual,yvisual,p{time_counter})
    %hold on
    %imagesc(xvisual,yvisual,p_desired)
    %plot(x,p_desired);

    xlim([-1 1]*xmaxvisual)
    ylim([-1 1]*xmaxvisual)
    %zlim([-1 1]*20e-6)
    %zlim([0 1]*0.05)


    xlabel('x')
    ylabel('y')
    zlabel('z')

    %hold off
    if VIDEO == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
    end

    if max(max(p{time_counter})) > 10
        break
    end


end

if VIDEO == 1
    close(VideoName);
end