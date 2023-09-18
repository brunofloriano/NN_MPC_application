clear all; close all; clc;

load('results/dataBalloon 2023-5-12-18-51/dataBalloon 2023-5-12-18-51.mat')
%load('results/dataBalloon 2023-5-13-6-17/dataBalloon 2023-5-13-6-17.mat')
%load('results/dataBalloon 2023-5-30-23-16/dataBalloon 2023-5-30-23-16.mat')

VideoName = plot_movie(savefolderlocal,VIDEO,GIF,t,individual_coord,Ggraph,GRAPH_AREA,Nr,encompassed_area_time);
