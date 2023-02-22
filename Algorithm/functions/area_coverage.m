function [area_covered, encompassed_area] = area_coverage(individual_coord,max_radius,map_steps, encompassed_area, area_radius,time_counter,total_circle_area)
N = length(individual_coord);
%radial_step = 1e3;
%angular_step = 2*pi/100;

plane_size = max_radius;
plane_step = 2*plane_size/map_steps;
%encompassed_area = zeros(map_steps, map_steps);

for agent_counter = 1:N
    x = individual_coord{agent_counter}(1,time_counter);
    y = individual_coord{agent_counter}(2,time_counter);

    % Check if new position is within bounds of the plane
    if x - area_radius < - plane_size
        x = -plane_size + area_radius;
    elseif x + area_radius > plane_size
        x = plane_size - area_radius;
    end
    if y - area_radius < - plane_size
        y = -plane_size + area_radius;
    elseif y + area_radius > plane_size
        y = plane_size - area_radius;
    end

    x_vector = (round((x-area_radius)/plane_step):round((x+area_radius)/plane_step)) + map_steps/2 + 1;
    y_vector = (round((y-area_radius)/plane_step):round((y+area_radius)/plane_step)) + map_steps/2 + 1;

    if sum(x_vector > map_steps) > 0 || sum(y_vector > map_steps) > 0
        encompassed_area(x_vector(1:end-1), y_vector(1:end-1)) = 1;
    else
        encompassed_area(x_vector, y_vector) = 1;
    end
  
end

encompassed_area = encompassed_area.*total_circle_area;

area_covered = sum(sum(encompassed_area))/sum(sum(total_circle_area)); %sum(sum(encompassed_area))*plane_step^2;

end