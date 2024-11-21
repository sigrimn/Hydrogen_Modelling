%% H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_PERTURBED_SCHEDULE_2

%% Parameter tuning of a very coarse upscaling of the H2 Benchmarking model  
% Optimering med matchObservedH2W_changeBC ser ut til å fungere. God
% optimering etter run2. Dette med de 3 øverste parameterne. VELDIG GODT
% RESULTAT!!!
% Når jeg kjører med alle parametre...:
clear all;
close all;
mrstModule add ad-core ad-blackoil deckformat agglom upscaling coarsegrid...
        mrst-gui ad-props incomp optimization test-suite linearsolvers 

%Define the name of the simulation for output files
name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT'; %changed name

%Read the Eclipse deck file containing the simulation data
deck = readEclipseDeck('/Users/sigridmarianese/Documents/MATLAB_workfolder/HydrogenOptimization/UHS/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');

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

% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
%  simulatePackedProblem(problem);

%%ALL above is for the original problem
name_perturbed = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_PERTURBED_SCHEDULE_2'; %changed name
problemRefPert = problem;
problemRefPert.BaseName = name_perturbed;

%% ONLY BHP CONTROL
ind = 1:100;
scheduleRefPert = problemRefPert.SimulatorSetup.schedule;
scheduleRefPert.step.val     = scheduleRefPert.step.val(ind);
scheduleRefPert.step.control = scheduleRefPert.step.control(ind);
    
dt       = scheduleRefPert.step.val;
nstep    = numel(scheduleRefPert.step.val);

% For the training run, we create a schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)

Wref    = scheduleRefPert.control(1).W;   %Should it here be all controls not only (1)?
scheduleRefPert_Bhp = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);

scheduleRefPert = scheduleRefPert_Bhp;
for i = 1:length(scheduleRefPert.control)
    scheduleRefPert.control(i).bc = schedule.control(i).bc;  
end


problemRefPert.SimulatorSetup.schedule = scheduleRefPert;

%% Define new output directory

% Define the new output directory
outputDir = '/Users/sigridmarianese/Documents/MATLAB_workfolder/Outputs/H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT_PERTURBED_SCHEDULE_2';

% Check and create the directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Create subdirectories for wellSols, states, and reports
mkdir(fullfile(outputDir, 'GenericBlackOilModel'));
mkdir(fullfile(outputDir, 'wellSols'));
mkdir(fullfile(outputDir, 'states'));
mkdir(fullfile(outputDir, 'reports'));
% here: simulatePackedProblem with new name and schedule
problemRefPert.OutputHandlers.wellSols.dataDirectory = outputDir;
problemRefPert.OutputHandlers.states.dataDirectory = outputDir;
problemRefPert.OutputHandlers.reports.dataDirectory = outputDir;

%% Simulate
simulatePackedProblem(problemRefPert) %, 'restartStep',1);

[wsRefPert, statesRefPert] = getPackedSimulatorOutput(problemRefPert);
%%

modelRefPert    = problemRefPert.SimulatorSetup.model;

%% Coarse-scale model
%Here we start building the coarse-scale model. This is to try to optimize
%the parametres for this model (but we want the boundary conditions for
%instance to match with the fine scaled model.
% We make a coarse grid defined by a uniform 6 x 6 x 1 partition 
blockIx = partitionUI(modelRefPert.G, [6,6,3]);
blockIx = processPartition(modelRefPert.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelCoarse = upscaleModelTPFA(modelRefPert, blockIx);
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
%figure, subplot(1,2,1)
%plotGrid(modelRef.G, 'EdgeAlpha',.2); 
%title('Fine-scale grid (18553 cells)')
%plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
%axis off tight, view(174,60), camlight headlight
%subplot(1,2,2)
%plotGrid(modelCoarse.G, 'EdgeAlpha',.8);
%plotWell(modelRef.G, Wref, 'Color', 'k', 'FontSize', 10);
%axis off tight, view(174,60), camlight headlight
%title('Coarse-scale grid (33 cells)')

% model.FacilityModel.ReservoirModel = wsRef;

% Initialize the reservoir solution based on the model grid and initial pressure
%state0 = initResSol(model.G, options.initPres, [1, 0]);        
%state0.rs = zeros(size(state0.pressure));  % Set initial solution gas-to-oil ratio to zero
% Initialize well solutions    
%wellSol = initWellSolAD(schedule.control(1).W, model, state0);
%wellSol.bhp = options.initPres;  % Set initial bottom hole pressure for all wells
%state0.wellSol = wellSol;  % Store well solutions in the state structure
%% Simulate initial upscaled coarse model and compare to reference

stateCoarse0   = upscaleState(modelCoarse, modelRefPert, state0); %rmfield(statesRef{index_initstate}, {'time', 'sMax', 'rv', 'status', 'FlowProps', 'PVTProps'}));  %initial state from reference at random time point 
scheduleCoarse = upscaleSchedule(modelCoarse, scheduleRefPert);
%%
modelCoarse = modelCoarse.validateModel(scheduleCoarse.control(1));
%%
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);
%%

plotWellSols({wsRefPert, wsCoarse}, ...
    {scheduleRefPert.step.val, scheduleCoarse.step.val},...
    'datasetnames',{'fine scale model','initial upscaled model'});

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
    'transmissibility', 1,        'log',                  [],     [1e-3,  1e3]
     'sgl',              1,     'linear',             [0, 1.0],              []
     'sgcr',             1,     'linear',             [0, 1.0],              []
     'sgu',              1,     'linear',             [.0, 1],              []
     'sogcr',            1,     'linear',             [0, 1.0],              []
     'krw',              1,     'linear',           [.0, 1.5],              []
     'kro',              1,     'linear',           [.0, 1.5],              []};
     
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

for i = 1:length(wsRefPert)
    bhp_ref=[bhp_ref; wsRefPert{i}.bhp];
    qGs_ref = [qGs_ref; wsRefPert{i}.qGs];
    qOs_ref = [qOs_ref; wsRefPert{i}.qOs];
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
wg = sum([scheduleRefPert.control(1).W.val]).*day;
qWs =[];
for i =1:length(wsRefPert)
 qWs(i) = abs(sum([wsRefPert{i}.qOs]));
end
ww = max(qWs).*day;
wp = (max(bhp_ref)-min(bhp_ref))/barsa;
% weighting  = {'WaterRateWeight',  1/(max(abs(qOs_ref))/day()), ...         %Gas (earlier water): hydrogen       #Change so that this is also dependent on number of time steps (different cycles)
%               'OilRateWeight',    1/(max(abs(qGs_ref))/day()),...     %%Oil: water 
%               'BHPWeight',        1/(max(abs(bhp_ref))*barsa)};

weighting  = {'WaterRateWeight',  1/(wg/day), ...
              'OilRateWeight',    1/(wg/day), ...
              'BHPWeight',        1/(wp*barsa)};
% make handle
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_changeBC(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRefPert);
% The calibration can be improved by taking a large number of iterations,

% but here we set a limit of 30 iterations
[v1, p_opt1, history1] = unitBoxBFGS(pvec, objh, 'objChangeTol', 1e-5, ...
                                  'maxIt', 10, 'logPlot', true);

%% Model calibration Levenberg-Marquard (using full Jacobian)
mismatchFn2 = @(model, states, schedule, states_ref, compDer, tstep, state) ...
   matchObservedH2W_changeBC(model, states, schedule, states_ref, ...
       'computePartials', compDer, 'tstep', tstep, weighting{:},...
       'state', state, 'from_states', false, 'mismatchSum', false);
objh2 = @(p) evaluateMatchSummands(p, mismatchFn2, setup_init, parameters, statesRefPert);
% The calibration can be improved by taking a large number of iterations,
% but here we set a limit of 30 iterations
[v2, p_opt2, history2] = unitBoxLM(p_opt1, objh2, 'maxIt', 10);

%% Saving the optimized variables

data.PERTURBED_SCHEDULE_2.run1.v1 = v1;
data.PERTURBED_SCHEDULE_2.run1.p_opt1 = p_opt1;
data.PERTURBED_SCHEDULE_2.run1.history1 = history1;
data.PERTURBED_SCHEDULE_2.run1.v2 = v2;
data.PERTURBED_SCHEDULE_2.run1.p_opt2 = p_opt2;
data.PERTURBED_SCHEDULE_2.run1.history2 = history2;
data.PERTURBED_SCHEDULE_2.run1.p_init = pvec;


% Save the structs 
save('PERTURBED_SCHEDULE_2_run1_paramsall.mat', '-struct', 'data');



%% Compare convergence history
%figure
semilogy(abs(history1.val(1:9))', '-o', 'LineWidth', 2);
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
%% compare reference, initial coarse and optimized coarse model outputs
plotWellSols({wsRefPert, wsCoarse, wellSols_opt1,wellSols_opt2}, ...
              repmat({scheduleRefPert.step.val}, 1, 4), ...
              'datasetnames', {'reference','coarse initial','coarse tuned (Q-N)','coarse tuned (L-M)'});

%% Plot the pore volume updates
% fetch pore volume differences in initial and tuned coarse models
oudpv = setup_opt1.model.operators.pv - setup_init.model.operators.pv;
figure
plotCellData(modelCoarse.G, oudpv, 'EdgeColor','none');
plotFaces(modelCoarse.G, boundaryFaces(modelCoarse.G), 'EdgeColor', [0.4 0.4 0.4], ...
         'EdgeAlpha',.5, 'FaceColor', 'none');
view(174,60);
plotWell(modelRefPert.G, Wref, 'Color', 'k', 'FontSize', 10); axis off tight
colorbar('south');
                        
%% 
% Compare reference, initial coarse and optimizes coarse for a different schedule 
rng(100);
W = schedule.control.W;



s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
for i = 1:length(s_new.control)
    s_new.control(i).bc = schedule.control(i).bc;
end
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRefPert, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt = setup_opt2.schedule.control(1).W;
bc_opt = setup_opt2.schedule.control(1).bc;
w_old = scheduleCoarse.control(1).W;
bc_old = scheduleCoarse.control(1).bc;   
s_opt = simpleSchedule(dt, 'W', w_opt,'bc',bc_opt);
s_old = simpleSchedule(dt, 'W', w_old,'bc',bc_old);
for kw = 1:numel(Wref)
    s_opt.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end
ws_coarse_new1 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt);
ws_coarse_old1 = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);

%%
plotWellSols({ws_new1,  ws_coarse_old1, ws_coarse_new1}, ...
    repmat({schedule.step.val(ind)}, 1, 3), ...
    'datasetnames', {'reference','coarse initial','coarse tuned'});

%% Final comparison
new_scheduleReference = scheduleReference;
for kw = 1:numel(scheduleReference.control)
    new_scheduleReference.control(kw).W(1).cstatus = w_opt.cstatus;
    new_scheduleReference.control(kw).W(1).fperf = w_opt.fperf;
    new_scheduleReference.control(kw).W(1).parentIndices = w_opt.parentIndices;
    new_scheduleReference.control(kw).W(1).r = w_opt.r;
    new_scheduleReference.control(kw).W(1).dir = w_opt.dir;
    new_scheduleReference.control(kw).W(1).rR = w_opt.rR;
    new_scheduleReference.control(kw).W(1).cells = w_opt.cells;
    new_scheduleReference.control(kw).W(1).cstatus = w_opt.cstatus;
    new_scheduleReference.control(kw).W(1).WI = w_opt.WI;
    new_scheduleReference.control(kw).W(1).dZ = w_opt.dZ;
    new_scheduleReference.control(kw).bc = bc_opt;
end
[wstuned1,statestuned1] = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, new_scheduleReference);

[wstuned2,statestuned2] = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, new_scheduleReference);


%% Running for entire time slot 2501 time steps. Compare with original reference model

%Define the name of the simulation for output files
name = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT'; %change this name to store another place  

%perturbated simulation

%Read the Eclipse deck file containing the simulation data
deck = readEclipseDeck('/Applications/mrst-2024a/modules/optimization/HydrogenOptimization/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');

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

% Pack the simulation problem with the initial state, model, and schedule
problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);
scheduleReference = problem.SimulatorSetup.schedule;

%Above copied from run_uhs_H2Brine_benchmark
%%
Wref     = schedule.control.W;
dt       = schedule.step.val;
nstep    = numel(schedule.step.val);

% Simulate
simulatePackedProblem(problem);

[wsRef, statesRef] = getPackedSimulatorOutput(problem);

%%
N = length(wsRef);
reference_bhp = zeros(N);
reference_qOs = zeros(N);
reference_qGs = zeros(N);

new_bhp = zeros(N);
new_qOs = zeros(N);
new_qGs = zeros(N);


for i = 1:2518
    reference_bhp(i) = wsRef{i}.bhp;
    reference_qOs(i) = wsRef{i}.qOs;
    reference_qGs(i) = wsRef{i}.qGs;

    new_bhp(i) = wstuned1{i}.bhp;
    new_qOs(i) = wstuned1{i}.qOs;
    new_qGs(i) = wstuned1{i}.qGs;
end
    
%%
% Plotting
%time = 1:N;  % Create a time vector corresponding to the length of A and B

% Plotting reference and fitted full models
figure; % Create a new figure
plot(reference_bhp, '-b');  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_bhp, '-r');  % Plot B with a red line
% Adding legends
legend('Reference model', 'Fitted model, L-M');

% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('barsa');  % Label for the y-axis
title('Comparison of reference model and fitteed model bottom hole pressure');  % Title for the plot



% Optional: Grid for better visualization
grid on;
hold off;
%%
% Plotting reference and fitted full models
figure; % Create a new figure
plot(reference_qOs, '-b', 'LineWidth', 1.5);  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_qOs, '-r', 'LineWidth', 1.5);  % Plot B with a red line

% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('m^3/day');  % Label for the y-axis
title('Comparison of reference model and fitteed model water flow');  % Title for the plot

% Adding legends
legend('Reference model', 'Fitted model, L-M');

% Optional: Grid for better visualization
grid on;
hold off;

%%
% Plotting reference and fitted full models
figure; % Create a new figure

plot(reference_qGs, '-b', 'LineWidth', 1.5);  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_qGs, '-r', 'LineWidth', 1.5);  % Plot B with a red line


% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('m^3/day');  % Label for the y-axis
title('Comparison of reference model and fitteed model hydrogen flow');  % Title for the plot

% Adding legends
legend('Reference model', 'Fitted model, L-M');


%legend('Reference model', 'Fitted model, L-M', 'b', 'r');
%ah1 = axes('position',get(gca,'position'),'visible','off');
%legend(ah1, [h1(2) h2(2)], {'Test1','Test2'}, 'Reference model', 'Fitted model, L-M');

% Optional: Grid for better visualization
grid on;
hold off;
%%
%plotWellSols({wsRef,  wstuned1, wstuned2}, ...
%    {new_scheduleReference.step.val, new_scheduleReference.step.val, new_scheduleReference.step.val},...
%    'datasetnames',{'reference model','fitted model 1', 'fitted model 2'});
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
