function area_4th = def_4th_area(Area)

x = Area.x/Area.delta_x;
y = Area.y/Area.delta_x;
k = 2;
area_4th = Area.encompassed_area;
sigma = 50e3/Area.delta_x;
variance = sigma^2;
covariance = variance*eye(2);
det_covariance = det(covariance);

for i = 1:Area.map_steps
    for j = 1:Area.map_steps
        vector = [x(i) y(j)]';
        %pdf = (2*pi)^(-k/2)*det_covariance^(-1/2)*exp( -1/2*norm(vector)^4);
        %pdf = (2*pi)^(-k/2)*det_covariance^(-1/2)*exp( -(x(i)^4 + y(j)^4 )/(2*sigma^4) );
        pdf = (2*pi)^(-k/2)*det_covariance^(-1/2)*exp( -( norm(vector)^4 )/(2*sigma^4) );
        area_4th(i,j) = pdf;
    end
end
max_point = max(max(area_4th));
area_4th = area_4th/max_point;

%mesh(area_4th)