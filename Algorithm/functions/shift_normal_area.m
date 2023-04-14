function normal_area = shift_normal_area(std_area,x_grid,y_grid,map_steps)

for i = 1:map_steps
    ishift = i + (round(map_steps/2) - x_grid);
    for j = 1:map_steps
        jshift = j + (round(map_steps/2) - y_grid);
        if ishift > map_steps || ishift < 1 || jshift > map_steps || jshift < 1
            normal_area(i,j) = 0;
        else
            normal_area(i,j) = std_area(ishift,jshift);
        end
    end
end