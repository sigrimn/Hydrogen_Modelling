%% H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_PERTURBED_SCHEDULE_2

%% Parameter tuning of a very coarse upscaling of the H2 Benchmarking model  
% Optimering med matchObservedH2W_changeBC ser ut til å fungere. God
% optimering etter run2. Dette med de 3 øverste parameterne. VELDIG GODT
% RESULTAT!!!
% Når jeg kjører med alle parametre...:

clear all;
close all;
addpath(genpath('/Applications/MRST/modules/H2store'));
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%Define the name of the simulation for output files

figureOutputFolder = '/Users/sigridmarianese/Documents/MATLAB_workfolder/project_figures/3D/Case2_new';
%Read the Eclipse deck file containing the simulation data

name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT'; %original model 
deck = readEclipseDeck('/Users/sigridmarianese/Documents/MATLAB_workfolder/HydrogenOptimization/UHS/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');


%% Prepare simulation parameters and initial state
problem = initialize3D(deck, name);
[wsRef_full, statesRef] = getPackedSimulatorOutput(problem);
%ALL above is for the original problem
name_perturbed = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_CASE2_NEW'; %changed name
problemRefPert = problem;
problemRefPert.BaseName = name_perturbed;

%% ONLY BHP CONTROL

scheduleRefPert = createPerturbedSchedule3D(problem.SimulatorSetup.schedule, 50, 50, 0, 0, .01, .25);

problemRefPert.SimulatorSetup.schedule = scheduleRefPert;

%% Define new output directory

% Define the new output directory
outputDir = '/Users/sigridmarianese/Documents/MATLAB_workfolder/Outputs/H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_CASE2_NEW';
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
      %'sgl',              1,     'linear',             [0, 1.0],              []
      %'sgcr',             1,     'linear',             [0, 1.0],              []
      %'sgu',              1,     'linear',             [.0, 1],              []
      %'sogcr',            1,     'linear',             [0, 1.0],              []
      %'krw',              1,     'linear',           [.0, 1.5],              []
      %'kro',              1,     'linear',           [.0, 1.5],              []};
     
parameters = setupParameters(setup_init, config);
%%
weighting = objectiveWeighting(wsRefPert);

% make handle
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_dynamic(model, states, schedule, states_ref,...
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
   matchObservedH2W_dynamic(model, states, schedule, states_ref, ...
       'computePartials', compDer, 'tstep', tstep, weighting{:},...
       'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRefPert);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 25 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt1, objh2, 'objChangeTol', 1e-5, 'maxIt', 25);

%% Save results
results = packResults(v1, p_opt1, history1, v2, p_opt2, history2, modelCoarse, modelRef, scheduleCoarse, scheduleRef, wsRef, wsCoarse, wellSols_opt1, wellSols_opt2);
save('3D_results_case2_3_BFGS_LM_dynamic_newtake.mat', 'results');

%% Run to plot the optimization
setup_opt1 = updateSetupFromScaledParameters(setup_init, parameters, p_opt1); 
[wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

% compare reference, initial coarse and optimized coarse model outputs
plotWellSols({wsRefPert, wsCoarse, wellSols_opt1,wellSols_opt2}, ...
              repmat({scheduleRefPert.step.val}, 1, 4), ...
              'datasetnames', {'reference','coarse initial','coarse tuned (BFGS)','coarse tuned (LM)'});

%% Re-run simulation for optimal parameters for full time-horizon

schedule = problem.SimulatorSetup.schedule;
n = 1400;
ind = 1:n;
scheduleshort.step.val     = schedule.step.val(ind);
scheduleshort.step.control = schedule.step.control(ind);
scheduleshort.control = schedule.control(1:schedule.step.control(n));

schedule = scheduleshort;

modelCoarse = makeCoarseModel3D(problem.SimulatorSetup.model, [6, 6, 3]);
[scheduleCoarse_init, stateCoarse0_init] = makeCoarseSchedule3D(modelCoarse, problem.SimulatorSetup.model,schedule, problem.SimulatorSetup.state0);

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

%% Coarse optimized
setup_opt = updateSetupFromScaledParameters(setup_init, parameters, p_opt2);
[scheduleCoarse, stateCoarse0] = makeCoarseSchedule3D(setup_opt.model, problem.SimulatorSetup.model, schedule, problem.SimulatorSetup.state0);
[ws_opt, states_opt] = simulateScheduleAD(setup_opt.state0, setup_opt.model, scheduleCoarse);


results.ws_init = ws_init;
results.ws_opt = ws_opt;
results.scheduleFull = schedule;

save('3D_results_case2_3_BFGS_LM_dynamic_newtake.mat', 'results');
%% Copyright Notice

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
