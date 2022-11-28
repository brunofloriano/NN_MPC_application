function output = generate_agents(N,max_radius)

matrix = rand(2,N);

angular_position = 2*pi*matrix(1,:);
radius = max_radius*matrix(2,:);

output_polar = [angular_position;radius];
[x,y] = pol2cart(angular_position,radius);
output_cart = [x;y];

output = output_cart;