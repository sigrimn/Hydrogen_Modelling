function result = optimize(problem, i, j)

simulatePackedProblem(problem, 'RestartStep', 1);
modelRef    = problem.SimulatorSetup.model;

[wsRef, statesRef] = getPackedSimulatorOutput(problem);

%% Create coarse model and schedule
modelCoarse = makeCoarseModel2D(modelRef, [11, 11]);  %not implemented change of dimensions
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule2D(modelCoarse, modelRef, problem.SimulatorSetup.schedule, problem.SimulatorSetup.state0);

[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

%% Set up misfit function and run optimization
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);

pv = modelCoarse.operators.pv;

config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.0001*pv, 10*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]      
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]};
     %'sgl',              1,     'linear',             [0, 1.0],              []
     %'sgcr',             1,     'linear',             [0, 1.0],              []
     %'sgu',              1,     'linear',             [.0, 1],              []
     %'sogcr',            1,     'linear',             [0, 1.0],              []
     %'krw',              1,     'linear',           [.0, 1.5],              []
     %'kro',              1,     'linear',           [.0, 1.5],              []};

parameters = setupParameters(setup_init, config);
weighting = objectiveWeighting(wsRef);


% make handle          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,

% but here we set a limit of 30 iterations
[v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
                                  'maxIt', 10, 'logPlot', false);

%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
   matchObservedH2W(model, states, schedule, states_ref, ...
       'computePartials', compDer, 'tstep', tstep, weighting{:},...
       'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt1, objh2, 'maxIt', 10, 'plotEvolution', false);

setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
[ws_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

result = struct();
result.v1 = v1;
result.v2 = v2;
result.p_opt1 = p_opt1;
result.p_opt2 = p_opt2;
result.history1 = history1;
result.history2 = history2;
result.i = i;
result.j = j;
result.wsRef = wsRef;
result.ws_opt2 = ws_opt2;
result.wsCoarse = wsCoarse;


