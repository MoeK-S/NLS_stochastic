%work out how to set default values in a function. Also should I really
%create a class and add the methods solve, and plot? that might be the best
%thing to do honestly. Shall we see how that would work at least?



solution = NLS_stochastic_solver;
solution.result = solve(solution);

%plot the sample path.
plotting(solution)