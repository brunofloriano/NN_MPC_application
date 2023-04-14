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
addpath functions/Bens_functions/

% Initialize time parameters
delta_t = 0.001;
tmax = 10;
t = 0:delta_t:tmax;

% Initialize space parameters
xmax = 3;
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

% p_initial = mvnpdf(xy,[0 0],[5 0;0 5]);
% p_initial = reshape(p_initial,length(x),length(x));


dist1.mean = [0;0]; dist1.std = [1;1]; dist1.P = [dist1.std(1,1)^2 0; 0 dist1.std(2,1)^2]; dist1.n = length(x)^2;
dist1.Pinv = inv(dist1.P); dist1.Pdet = det(dist1.P);
% 
% dist2.mean = [-3;-3]; dist2.std = [0.2;0.2]; dist2.P = [dist2.std(1,1)^2 0; 0 dist2.std(2,1)^2]; dist2.n = round(((6*dist2.std(1,1)/G.dx)+1)*((6*dist2.std(2,1)/G.dx)+1));
% dist2.Pinv = inv(dist2.P); dist2.Pdet = det(dist2.P); 

%Hole in the middle Distribution
D_t = Initialize_PDF(G,dist1); D_t.P = D_t.P-D_t.P(2);

D = my_initial_PDF(p_initial,G); D = Initialize_vuw(D,G);

%D = Initialize_PDF(G,dist1); D = Initialize_vuw(D,G); 
% D_ag = Initialize_PDF(G,dist2); D = combine_PDF(D,D_ag,G);
% 
% for i=1:(2*3/G.dx)
%     dist2.mean = dist2.mean + [G.dx; G.dx]; 
%     D_ag = Initialize_PDF(G,dist2); 
%     D = combine_PDF(D,D_ag,G); 
% end

D = boundary_conditions(D,G); plot_initial(D,G);
figure
K = RHS_P(D,G); 
% D_new.P = dictionary(); 
% D_new.P=D.P+G.dt.*K;
% diff = sum(abs(D.P - D_new.P)); 
diff = 1e10;
timestep = 1; clear F_top; clear F_dtv; %clear F_side;  clear F_iso; 
while(diff > 0.02)
  k1=RHS_P(D,G); D.P = D.P + (G.dt/2).*k1;
  k2=RHS_P(D,G); D.P = D.P + (G.dt/2).*k2;
  k3=RHS_P(D,G); D.P = D.P + G.dt.*k3;
  k4=RHS_P(D,G);
  D_new.P=D.P+(G.dt/6).*k1+(G.dt/3).*(k2+k3)+(G.dt/6).*k4;
  %[diff, F1, F2] = plot_PDF(D,G,timestep,dist1,D_t,D_new);
  diff = sum(abs(D.P - D_new.P)); 
  %F_top(timestep) = F1; F_dtv(timestep) = F2;
  p = D2p_convert(D);
  mesh(x,y,p)
  timestep=timestep+1; D.P = D_new.P;  

      %hold off
    if VIDEO == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
    end

    if max(max(p)) > 10
        break
    end

end

if VIDEO == 1
    close(VideoName);
end