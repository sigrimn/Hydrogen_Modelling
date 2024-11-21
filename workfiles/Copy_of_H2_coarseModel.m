%% Parameter tuning of a very coarse upscaling of the H2 Benchmarking model   
% Kun en coarse modell. 100 time steps. 
% Med matchObserved_changeBC som er konstant weighting, så oppnår denne
% IKKE konvergens etter første iterasjon, men etter andre iterasjon når man
% 1e-3-objective-verdi (nesten 10^-4). Lagret figur: Coarse_model_static_weights_params123.fig'
% Med matchObserved 
close all;
clear all;
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%Define the name of the simulation for output files
name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_COARSE'; %change this name to store another place  

%Read the Eclipse deck file containing the simulation data
deck = readEclipseDeck('/Users/sigridmarianese/Documents/MATLAB_workfolder/HydrogenOptimization/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');
%% Prepare simulation parameters and initial state
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

problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

% We focus on the first 100 steps
ind = 1:100;
schedule = problem.SimulatorSetup.schedule;
schedule.step.val     = schedule.step.val(ind);
schedule.step.control = schedule.step.control(ind);
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

schedule.control = schedule.control(1:schedule.step.control(100));

% For the training run, we create a schedule with varying controls
%perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
%rng(0)
%schedule = perturbedSimpleSchedule(dt, 'W', Wref, ...
%    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
% Pack the simulation problem with the initial state, model, and schedule

% Simulate

problem.SimulatorSetup.schedule = schedule;
simulatePackedProblem(problem);

[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;

%Here we start building the coarse-scale model. This is to try to optimize
%the parametres for this model (but we want the boundary conditions for
%instance to match with the fine scaled model.
%% Coarse-scale model
% We make a coarse grid defined by a uniform 6 x 6 x 1 partition 
x = 6;
y = 6;
z = 3;
blockIx = partitionUI(modelRef.G, [x, y, z]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).
pts = modelCoarse.fluid.krPts;
scaling = {'SGL',   pts.g(1), 'SGCR', pts.g(2), 'SGU', pts.g(3), ...
            'SOGCR', pts.og(2), 'KRG',  pts.g(4), 'KRO', pts.og(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

%% Plot reference and coarse model grids with wells
figure, subplot(1,2,1)
plotGrid(modelRef.G, 'EdgeAlpha',.2); 
title('Fine-scale grid (18553 cells)')
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
subplot(1,2,2)
plotGrid(modelCoarse.G, 'EdgeAlpha',.8);
plotWell(modelRef.G,  Wref, 'Color', 'k', 'FontSize', 10);
axis off tight, view(174,60), camlight headlight
title(sprintf('Coarse-scale grid %d cells', x*y*z));

%% Simulate initial upscaled coarse model and compare to reference
stateCoarse0   = upscaleState(modelCoarse, modelRef, state0);
scheduleCoarse = upscaleSchedule(modelCoarse, schedule, 'wellUpscaleMethod', 'sum');
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

plotWellSols({wsRef, wsCoarse}, ...
    {schedule.step.val, scheduleCoarse.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'});

%% Specify parameters for tuning
setup_init = struct('model', modelCoarse, 'schedule', scheduleCoarse, 'state0', stateCoarse0);
pv = modelCoarse.operators.pv;
% set up 'matrix' for parameter options for easier editing. The specific
% limits set for the various parameters influences the tuning/optimization 
% procedure to a large extent
config = {...
     %name           include    scaling              boxlims   relativeLimits  
    'porevolume',       1,     'linear',    [.001*pv, 3*pv],              []   
    'conntrans',        1,        'log',                  [],     [1e-3, 1e3]      
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]  
    'sgl',              1,     'linear',             [0, .3],              []
    'sgcr',             1,     'linear',             [0, .4],              []
    'sgu',              1,     'linear',             [.7, 1],              []
    'sogcr',            1,     'linear',             [0, .4],              []
    'krw',              1,     'linear',           [.5, 1.5],              []
    'kro',              1,     'linear',           [.5, 1.5],              []};
    
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    parameters = addParameter(parameters, setup_init, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits',config{k,5});
end
%%
bhp_ref = [];
qGs_ref = [];
qOs_ref = []; 

for i = 1:length(wsRef)
    bhp_ref=[bhp_ref; wsRef{i}.bhp];
    qGs_ref = [qGs_ref; wsRef{i}.qGs];
    qOs_ref = [qOs_ref; wsRef{i}.qOs];
end 

%% Define the mismatch function
% Function weighting influences the match of each quantity. Rate-weighting
% should be on the same order as (inverse of) rates. BHP-weighting on the
% order of pressure drop in the model.
%bhp_relative_weight = (max(bhp_ref) - max(bhp))/max(bhp_ref);
%qGs_relative_weight = (max(qGs_ref) - max(qGs))/max(qGs_ref);
%qOs_relative_weight = (max(abs(qOs_ref)) - max(abs(qOs)))/max(abs(qOs_ref));
%weighting  = {'WaterRateWeight',  1/(qGs_relative_weight/day()), ...         %Gas (earlier water): hydrogen       #Change so that this is also dependent on number of time steps (different cycles)
%              'OilRateWeight',    1/qGs_relative_weight,...     %/(qOs_relative_weight/day), ... ...   %Oil: water 
%              'BHPWeight',        1/(bhp_relative_weight*barsa)};

weighting  = {'WaterRateWeight',  1/(max(abs(qOs_ref))/day()), ...         %Gas (earlier water): hydrogen       #Change so that this is also dependent on number of time steps (different cycles)
              'OilRateWeight',    1/(max(abs(qGs_ref))/day()),...     %%Oil: water 
              'BHPWeight',        1/(max(abs(bhp_ref))*barsa)};


% make handle          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_changeBC(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef(ind));
% The calibration can be improved by taking a large number of iterations,

% but here we set a limit of 30 iterations
[v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
                                  'maxIt', 10, 'logPlot', true);

%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_changeBC(model, states, schedule, states_ref, ...
        'computePartials', compDer, 'tstep', tstep, weighting{:},...
        'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRef);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt1, objh2, 'maxIt', 10);

%%
% Store the values in a file
% Save variables to 'data.mat'
% Create a struct with custom field names
data_121124.CoarseModel_allparams_121124.Run1.v2 = v2;              %run2 are with all params
data_121124.CoarseModel_allparams_121124.Run1.p_opt2 = p_opt2;
data_121124.CoarseModel_allparams_121124.Run1.history2 = history2;

% Save the structs 
save('CoarseModel_tuned_on_new_schedule_allparams_Run1_121124.mat', '-struct', 'CoarseModel_allparams_121124');


%% Compare convergence history
figure
semilogy(abs(history1.val)', '-o', 'LineWidth', 2);
hold on
semilogy(abs(history2.val)', '--o', 'LineWidth', 2);
set(gca, 'Fontsize', 12); grid on
xlabel('Iteration'), ylabel('Residual')
legend({'Quasi-Newton'}) %, 'Levenberg-Marquard'})

%% Create new coarse model setups with the optimized parameters, and rerun 
%  for the optimized parameters
setup_opt1 = updateSetupFromScaledParameters(setup_init, parameters, p_opt1); 
[wellSols_opt1, states_opt1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, setup_opt1.schedule);
setup_opt2 = updateSetupFromScaledParameters(setup_init, parameters, p_opt2); 
[wellSols_opt2, states_opt2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);
%compare reference, initial coarse and optimized coarse model outputs

%%
plotWellSols({wsRef, wsCoarse, wellSols_opt2}, ...
              repmat({schedule.step.val}, 1, 3), ...
              'datasetnames', {'reference','coarse initial','coarse tuned (L-M)'});

%% Plot the pore volume updates
% fetch pore volume differences in initial and tuned coarse models
oudpv = setup_opt2.model.operators.pv - setup_init.model.operators.pv;
figure
plotCellData(modelCoarse.G, oudpv, 'EdgeColor','none');
plotFaces(modelCoarse.G, boundaryFaces(modelCoarse.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60);
plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
colorbar('south');
                        
%% Compare reference, initial coarse and optimizes coarse for a different schedule 
rng(100);
W = schedule.control.W;

dt = dt(1:10);

s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt = setup_opt2.schedule.control(1).W;
w_old = scheduleCoarse.control(1).W;   
s_opt = simpleSchedule(dt, 'W', w_opt);
s_old = simpleSchedule(dt, 'W', w_old);
for kw = 1:numel(Wref)
    s_opt.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end
ws_coarse_new = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt);
ws_coarse_old = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);

%%

plotWellSols({ws_new1,  ws_coarse_old, ws_coarse_new}, ...
    repmat({dt}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'});

savefig('Coarse_model_static_weights_params123.fig');

%% Here compare reference model, initial coarse and tuned coarse. Maybe for full time schedule? Maybe downscale?

%[w_new_coarse, states_new_coarse] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, setup_opt2.schedule);

%% Copyright Notice
%
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
