clear; close all; clc;

addpath functions
addpath functions/Bens_functions/

addpath 'Balloon'

start
save_folder = 'results';

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


plane_size = max_radius;
plane_step = 2*plane_size/map_steps;
%area_covariance_matrix = (25e3)^2*eye(2);
covariance = round(area_covariance_matrix/plane_step^2);
G.dt = 0.001;

normal_area = def_normal_area(map_steps,[30;30] + [1;1]*0,covariance);
for i = 1:10:1000
new_normal_area = def_normal_area(map_steps,[30;30] + [1;1]*i,covariance);
sum_area = normal_area;
sum_area(normal_area < new_normal_area) = new_normal_area(normal_area < new_normal_area);
normal_area = sum_area;
%prob_sum = (G.dx)^G.d.*sum(sum(normal_area)); 
%normal_area = normal_area/prob_sum;% D.n = length(D.P);
relax_area = relaxing_PDF(normal_area,G);


mesh(x,y,normal_area)
%imagesc(normal_area)
%mesh(x,y,relax_area)
%imagesc(relax_area)

 if VIDEO == 1
        VideoFrame = getframe(gcf);
        writeVideo(VideoName,VideoFrame);
    end


normal_area = relax_area;

end


if VIDEO == 1
    close(VideoName);
end 