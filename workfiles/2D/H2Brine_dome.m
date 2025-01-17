clear all;
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr spe10 coarsegrid
addpath(genpath('/Applications/MRST/modules/H2store'))
addpath(genpath('/Users/sigridmarianese/Documents/ProsjektOppgave'))
% Denne ble kjørt lørdag kveld, og fungerte veldig godt. Perturbering på
% injection og bhp control
%% Define the case name and read the Eclipse deck file
name = 'H2_STORAGE_RS_ILLUSTRATIVE_2D';

%Change name to fit exp
name_perturbed = 'H2_STORAGE_RS_COARSEMODEL_EXPX';
%% Use H2STORAGE_RS_SALT.DATA for brine
deck = readEclipseDeck('/Applications/MRST/modules/H2store/data/Illustrative_example/H2STORAGE_RS.DATA');

%% Set up the simulation parameters and model components
[~, ~, state0, model, schedule, ~] = H2Brine_initialization(deck);

%% Set up the linear and nonlinear solvers

lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                    % Assign the linear solver to the nonlinear solver

problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
%% Create perturbed reference schedule
% Change schedule to fit exp
scheduleReference = createPerturbedSchedule2D(problem.SimulatorSetup.schedule, 0, 100, 0, 0, 0, 0.4);  % Tror at man ikke måååå ha bhp her

%% Simulate Reference schedule
problem.SimulatorSetup.schedule = scheduleReference;
simulatePackedProblem(problem);
modelRef    = problem.SimulatorSetup.model;

[wsRef, statesRef] = getPackedSimulatorOutput(problem);

%% Create coarse model and schedule
modelCoarse = makeCoarseModel2D(modelRef, [11, 11]);  %not implemented change of dimensions
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule2D(modelCoarse, modelRef, scheduleReference, state0);

%% Plot
%plotWell2d(modelRef, modelCoarse, scheduleReference, scheduleCoarse);

%% Simulate coarse schedule
modelCoarse = modelCoarse.validateModel();
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

plotWellSols({wsRef, wsCoarse}, ...
    {scheduleReference.step.val, scheduleCoarse.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'});

%% Specify parameters for tuning
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
pv = modelCoarse.operators.pv;

config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.0001*pv, 10*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]      
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]
     'sgl',              1,     'linear',             [0, 1.0],              []
     'sgcr',             1,     'linear',             [0, 1.0],              []
     'sgu',              1,     'linear',             [.0, 1],              []
     'sogcr',            1,     'linear',             [0, 1.0],              []
     'krw',              1,     'linear',           [.0, 1.5],              []
     'kro',              1,     'linear',           [.0, 1.5],              []};

%parameters = setupParameters(setup_init, config);
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, setup_init, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits',config{k,5});
end

%% Set up misfit function and run optimization
% weights
weighting = objectiveWeighting(wsRef);


% make handle          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_dynamic(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

% %% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
%% The calibration can be improved by taking a large number of iterations,
% 
% % but here we set a limit of 30 iterations
[v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 4e-5, ...
                                  'maxIt', 10, 'logPlot', true);

%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
   matchObservedH2W_static(model, states, schedule, states_ref, ...
       'computePartials', compDer, 'tstep', tstep, weighting{:},...
       'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(pvec, objh2, 'maxIt', 30);

%% Create new coarse model setups with the optimized parameters, and rerun 
setup_opt1 = updateSetupFromScaledParameters(setup_init, parameters, p_opt1); 
[wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

% Plot
plotWellSols({wsRef, wsCoarse, wellSols_opt1, wellSols_opt2}, ...
              repmat({scheduleReference.step.val}, 1, 4), ...
              'datasetnames', {'reference','coarse initial', 'coarse tuned (BFGS)', 'coarse tuned (LM)'});

%% Pack and save results
results = packResults(v1, p_opt1, history1, v2, p_opt2, history2, modelCoarse, modelRef, scheduleCoarse, scheduleReference);

resultsFolder = "/Users/sigridmarianese/Documents/ProsjektOppgave/2D_results_case2/";
save(resultsFolder + "EXP1.mat", 'results');

%% Re-run simulation for optimal parameters for full time-horizon
[~, ~, state0, model, schedule, ~] = H2Brine_initialization(deck);
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
[wsRef_full, statesRef] = getPackedSimulatorOutput(problem);

%chose length of simulation
n = 1400;
ind = 1:n;
scheduleshort.step.val     = schedule.step.val(ind);
scheduleshort.step.control = schedule.step.control(ind);
scheduleshort.control = schedule.control(1:schedule.step.control(n));

schedule = scheduleshort;

%% Coarse initial
[scheduleCoarse_init, stateCoarse0_init] = makeCoarseSchedule2D(setup_init.model, model, schedule, state0);
[ws_init, states_init] = simulateScheduleAD(setup_init.state0, setup_init.model, scheduleCoarse_init);

%% Coarse optimized only BFGS
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, results.p_opt1);
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule2D(setup_opt.model, model, schedule, state0);
[ws_opt_BFGS, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, scheduleCoarse);

%% Coarse optimized LM
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, results.p_opt2);
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule2D(setup_opt.model, model, schedule, state0);
[ws_optLM, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, scheduleCoarse);
%% 
results.wsRef_full = wsRef_full;
results.ws_opt_BFGS = ws_opt_BFGS;
results.ws_opt_LM = ws_optLM;
results.ws_init = ws_init;
save(resultsFolder + "EXP1.mat", 'results');