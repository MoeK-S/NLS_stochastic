%work out how to set default values in a function. Also should I really
%create a class and add the methods solve, and plot? that might be the best
%thing to do honestly. Shall we see how that would work at least?


h = 0.00002;
no_timesteps = 1000000;
delta = 0.0002;

%the optional last 2 arguments can be used to specify a seed.
solution = NLS_stochastic_solver(h, no_timesteps, delta, 13, 24);
solution.result = solve(solution);

%plot the sample path.
plotting(solution)