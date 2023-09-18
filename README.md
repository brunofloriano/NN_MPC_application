# NN_MPC
Coordinating a Hurricane Monitoring System Using a Swarm of Buoyancy-Controlled Balloons Trading Off Communication and Coverage by using Neural network-based Model Predictive Control (NNMPC).

By Bruno Floriano, Benjamin Hanson, Thomas Bewley, Joao Ishihara and Henrique Ferreira (2023).

A link to the paper will be added here when published.

## How to run algorithm
The main file is "main.m". There are 3 scenarios that were simulated in the paper:

* Simulation 1 with T=6h, N=9 (default)
* Simulation 2 with T=12h, N=9
* Simulation 3 with T=6h, N=4

The default simulation is Simulation 1. To run Simuation 2, uncomment line 8 of the file 'functions/accumulate_laplacian.m'. To run Simulation 3, uncomment line 6 of the file 'Balloon/start.m'.

The output results in 2 videos of the system's evolution over time:

* Video 1: Top view of the balloons in a colormap of the interest function;
* Video 2: 3D plot of the interest function.

## Results

The folder "results" contain all the data (.mat) shown in the Results section. It is organized according to the folder's names as follows:

* Simulation 1: dataBalloon 2023-5-12-18-51;
* Simulation 2: dataBalloon 2023-5-13-6-17;
* Simulation 3: dataBalloon 2023-5-30-23-16.

### Videos

To generate the frames and videos of each simulation, run "make_movie.m". Uncomment the following lines for each simulation:

* Line 3 to obtain the videos of the Simulation 1;
* Line 4 to obtain the videos of the Simulation 2;
* Line 5 to obtain the videos of the Simulation 3.

### Graphs

To generate the area coverage/communication energy graph, run "analysis_balloons.m" (default is Simulation 1). For Simulation 2, uncomment line 16.
To make Simulation 3, and see the laplacian entries over time, uncomment line 17.
