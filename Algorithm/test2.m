% Clear workspace
%clear; clc; close all;
addpath Balloon/
addpath functions
start

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
addpath functions/Bens_functions/


% Initialize time parameters
delta_t = 0.001;
tmax = 10;
t = 0:delta_t:tmax;

% Initialize space parameters
xmax = 10;
%delta_x = 0.1; x_steps = 2*xmax/delta_x + 1;
x_steps = 61; delta_x = 2*xmax/(x_steps - 1);
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
D_FP = Dn*[1 0;0 1];
%D = Dn*[1 1;1 1];

G.dx=delta_x; G.dt=.001; G.xh=G.dx/2; G.d=2; G.dif = 1; 
G.x = x; G.y = y;

%% Initialize PDF as a uniform distribution

p_initial = zeros(n);
inside = p_initial;
round_boundary = inside;
un = p_initial;
uniform_lim = 2;
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
        un(xcounter,ycounter,1:2) = -D_FP*r^3*[cos(theta); sin(theta)];
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

% p_initial = mvnpdf(xy,[0 0],[0.1 0;0 0.1]);
% p_initial = reshape(p_initial,length(x),length(x));

 p_initial = def_4th_area(Area);
%p_initial = p_desired;

p = p_initial; diff = 1e10;
%D = p2D_convert(p,G);

while(diff > 0.002)
    mesh(x,y,p)
    drawnow
    %D_new = relaxing_PDF(D,G);
    p_new = test_relaxing_PDF(p);
    %diff = max(max(abs(p - p_new)));
    %D = D_new;
    %p = D2p_convert(D);
    p = p_new/max(max(p_new));

      %hold off
    if VIDEO == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
    end


end

if VIDEO == 1
    close(VideoName);
end