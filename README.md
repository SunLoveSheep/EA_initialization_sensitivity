# EA_initialization_sensitivity

Here are .cpp and .h files for simulation program that investigate EA's sensitivity to initial solutions.

Core idea is that, with a set of initial solution, EA can run to reach an optimization result. If we know the optima, how about move 1 
initial solution into the neighborhood of the optima? How will this affect the EA performance? How about move 2?

By adjusting number of initial solution moved and the range of "neighborhood", we try to collect data by this program and plot some curves
to investigate the influences from initial solutions.
