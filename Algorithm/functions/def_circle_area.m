%% Define circle area
% This function creates the matrix 'total_circle_area' with size
% 'map_steps'. It represents a map with '1' where it is inside a circle of 
% diameter 'map_steps' and '0' otherwise.

function total_circle_area = def_circle_area(map_steps)

total_circle_area = ones(map_steps);
radius = map_steps/2;

for i = 1:map_steps
    real_x = i - 1 - round(radius);
    for j = 1:map_steps
        real_y = j - 1 - round(radius);
        norm_radius = norm([real_x; real_y]);
        if norm_radius > radius
            total_circle_area(i,j) = 0;
        end
    end
end