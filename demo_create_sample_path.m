%this file demonstrates the use of the numerical algorithm to generate one
%sample path of the stochastic process X which solves the SDE given by
%Equation (7.1) in https://arxiv.org/abs/2205.12868 

h = 0.00002;
no_timesteps = 1000000;
delta = 0.0002;

%the optional last 2 arguments can be used to specify a seed.
solution = NLS_stochastic_solver(h, no_timesteps, delta);
solution.result = solve(solution);

%plot the sample path.
plotting(solution)