%% Define normal area
% This function creates the matrix 'normal_area' with size
% 'map_steps'. It represents a map with gaussian distribution centered at
% 'mean' and with covariance matrix 'covariance'.

function normal_area = def_normal_area(map_steps,mean,covariance)

%normal_area = zeros(map_steps);
%k = 2;
%det_covariance = det(covariance);
%inv_covariance = inv(covariance);
[x,y] = meshgrid(1:map_steps,1:map_steps);
X = [x(:) y(:)];
result = mvnpdf(X,mean',covariance);
normal_area = reshape(result,length(y),length(x));
max_area = max(max(normal_area));
normal_area = normal_area/max_area;

% for x = 1:map_steps
%     x_shift = x - mean(1);
%     for y = 1:map_steps
%         y_shift = y - mean(2);
%         vector = [x_shift;y_shift];
% 
%         pdf = (2*pi)^(-k/2)*det_covariance^(-1/2)*exp( -1/2*vector'*inv_covariance*vector);
% 
%         normal_area(x,y) = pdf;
% 
%     end
% end


% radius = map_steps/2;
% 
% for i = 1:map_steps
%     real_x = i - 1 - round(radius);
%     for j = 1:map_steps
%         real_y = j - 1 - round(radius);
%         norm_radius = norm([real_x; real_y]);
%         if norm_radius > radius
%             normal_area(i,j) = 0;
%         end
%     end
% end