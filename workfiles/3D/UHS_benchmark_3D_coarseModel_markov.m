%% H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_PERTURBED_SCHEDULE_2
%% Before running this, make a results folder in the markov work directory
%% start tmux session by writing desired name_perturbed (with double fnutts "  not '), rateFac and pressureFac. ex:
% name_perturbed = "UHS_005_001
% rateFac = 0.05
% pressureFac = 0.01

addpath(genpath('/work/sigrimn/MRST/modules/H2store'));

mrstModule add ad-core
mrstModule add ad-blackoil
mrstModule add deckformat
mrstModule add agglom
mrstModule add upscaling
mrstModule add coarsegrid
mrstModule add mrst-gui
mrstModule add ad-props
mrstModule add incomp
mrstModule add optimization
mrstModule add test-suite
mrstModule add linearsolvers



name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT'; %original model 
deck = readEclipseDeck('/work/sigrimn/MRST/modules/H2store/data/uhs_benchmark/UHS_BENCHMARK_RS.DATA');


%% Prepare simulation parameters and initial state
problem = initialize3D(deck, name);

%% ONLY BHP CONTROL

scheduleRefPert = createPerturbedSchedule3D(problem.SimulatorSetup.schedule, 50, 50, 0, 0, pressureFac, rateFac);  %Changed so that input is from terminal

problemRefPert.SimulatorSetup.schedule = scheduleRefPert;

%% Define new output directory
% Define the new output directory
outputDir = ["/work/sigrimn/Outputs/", name_perturbed]; 
problemRefPert = createOutputDirectory(outputDir, problemRefPert);


%% Simulate
simulatePackedProblem(problemRefPert);
[wsRefPert, statesRefPert] = getPackedSimulatorOutput(problemRefPert);
%% 
modelRefPert    = problemRefPert.SimulatorSetup.model;
%% Coarse-scale model
modelCoarse = makeCoarseModel3D(modelRefPert, [6, 6, 3]);
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule3D(modelCoarse, modelRefPert, scheduleRefPert, problemRefPert.SimulatorSetup.state0);

%% Simulate initial upscaled coarse model and compare to reference
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

%% Specify parameters for tuning
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
pv = modelCoarse.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent
config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.0001*pv, 10*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]      
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]};
     % 'sgl',              1,     'linear',             [0, 1.0],              []
     % 'sgcr',             1,     'linear',             [0, 1.0],              []
     % 'sgu',              1,     'linear',             [.0, 1],              []
     % 'sogcr',            1,     'linear',             [0, 1.0],              []
     % 'krw',              1,     'linear',           [.0, 1.5],              []
     % 'kro',              1,     'linear',           [.0, 1.5],              []};
     
parameters = setupParameters(setup_init, config);
%%
weighting = objectiveWeighting(wsRefPert);

% make handle
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRefPert);
% The calibration can be improved by taking a large number of iterations,

% but here we set a limit of 25 iterations
[v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 5e-4, ...
                                  'maxIt', 15, 'updateTol', 1e-5, 'logPlot', true);

%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
   matchObservedH2W(model, states, schedule, states_ref, ...
       'computePartials', compDer, 'tstep', tstep, weighting{:},...
       'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRefPert);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 25 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt1, objh2, 'objChangeTol', 5e-4, 'maxIt', 20);



%% Run to plot the optimization

setup_opt1 = updateSetupFromScaledParameters(setup_init, parameters, p_opt1); 
[wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

%% Save results
results = packResults(v1, p_opt1, history1, v2, p_opt2, history2, modelCoarse, modelRef, scheduleCoarse, scheduleRef, wsRefPert, wsCoarse, wellSols_opt1, wellSols_opt2);

%% Re-run simulation for optimal parameters for full time-horizon
% Laster variablene jeg trenger på nytt:
[v1, p_opt1, history1, v2, p_opt2, history2, modelCoarse, modelRef, scheduleRef, scheduleCoarse] = loadResults(results);

name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT'; %original model 
deck = readEclipseDeck('/work/sigrimn/MRST/modules/H2store/data/uhs_benchmark/UHS_BENCHMARK_RS.DATA');


%% Prepare simulation parameters and initial state
problem = initialize3D(deck, name);
[wsRef_full, statesRef] = getPackedSimulatorOutput(problem);

%%
modelCoarse = makeCoarseModel3D(problem.SimulatorSetup.model, [6, 6, 3]);
[scheduleCoarse_init, stateCoarse0_init] = makeCoarseSchedule3D(modelCoarse, problem.SimulatorSetup.model, problem.SimulatorSetup.schedule, problem.SimulatorSetup.state0);

setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse_init, 'state0', stateCoarse0_init);
[ws_init, states_init] = simulateScheduleAD(setup_init.state0, setup_init.model, scheduleCoarse_init);


setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse_init, 'state0', stateCoarse0_init);
pv = modelCoarse.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent
config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.0001*pv, 10*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]      
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]};
     % 'sgl',              1,     'linear',             [0, 1.0],              []
     % 'sgcr',             1,     'linear',             [0, 1.0],              []
     % 'sgu',              1,     'linear',             [.0, 1],              []
     % 'sogcr',            1,     'linear',             [0, 1.0],              []
     % 'krw',              1,     'linear',           [.0, 1.5],              []
     % 'kro',              1,     'linear',           [.0, 1.5],              []};
     
parameters = setupParameters(setup_init, config);

%Coarse optimized
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt2);
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule3D(setup_opt.model, problem.SimulatorSetup.model, problem.SimulatorSetup.schedule, problem.SimulatorSetup.state0);
[ws_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, scheduleCoarse);

results.ws_init = ws_init;
results.ws_opt = ws_opt;
save(["/work/sigrimn/results/", name_perturbed, ".mat"], 'results');

%% 

% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
