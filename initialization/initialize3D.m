function problem = initialize3D(deck, name)


[~, options, state0, model, schedule, ~] = modified_uhs_benchmark(deck);

% Add custom output functions to the model for additional diagnostics
model.OutputStateFunctions{end + 1} = 'CapillaryPressure';  % Output capillary pressure
model.OutputStateFunctions{end + 1} = 'SurfaceDensity';     % Output surface density
model.OutputStateFunctions{end + 1} = 'ShrinkageFactors';    % Output shrinkage factors
model.outputFluxes = false;                                   % Disable output of fluxes

%% Set up the linear and nonlinear solvers
lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                    % Assign the linear solver to the nonlinear solver

% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);