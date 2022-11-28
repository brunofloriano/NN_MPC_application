clear; close all; clc;
addpath 'functions'

N = 100;
max_radius = 200e3;
comm_range = 50e3;
max_root_connections = 5;

agents_position = generate_agents(N,max_radius);
position = [[0;0] agents_position];

%Ggraph = graph(ones(1,N),2:N+1);
Ggraph = routing_protocol(position,comm_range,max_root_connections);

plot(Ggraph,XData=position(1,:),YData=position(2,:))

%scatter(position(1,:),position(2,:))
