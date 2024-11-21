clearvars; 
mrstModule add ad-core ad-blackoil ad-props deckformat mrst-gui upr spe10 coarsegrid

%% Define the case name and read the Eclipse deck file
name = 'H2_STORAGE_RS_COARSEMODEL_PERTURBEDSCHEDULE_LOW_TOL';
%% Use H2STORAGE_RS_SALT.DATA for brine
deck = readEclipseDeck('/Applications/MRST/modules/H2store/data/Illustrative_example/H2STORAGE_RS.DATA');

%% Set up the simulation parameters and model components
[~, ~, state0, model, schedule, ~] = H2_illustration_storage_example_cartesian_grid(deck);
state0=rmfield(state0,'wellSol');
%%
% Add custom output functions to the model for additional diagnostics
model.OutputStateFunctions{end + 1} = 'CapillaryPressure';  % Output capillary pressure
model.OutputStateFunctions{end + 1} = 'SurfaceDensity';     % Output surface density
model.OutputStateFunctions{end + 1} = 'ShrinkageFactors';    % Output shrinkage factors
model.outputFluxes = false;                                   % Disable output of fluxes
model.toleranceCNV = 1.0e-5;
% model.FacilityModel.toleranceWellBHP = 1000;
% model.FacilityModel.toleranceWellRate = 1e-7;
% model.FacilityModel.toleranceWellMS = 1e-7;
%% Set up the linear and nonlinear solvers

lsolve = selectLinearSolverAD(model);                          % Select the linear solver for the model
nls = NonLinearSolver();                                       % Create a nonlinear solver object
nls.LinearSolver = lsolve;                                    % Assign the linear solver to the nonlinear solver

problem = packSimulationProblem(state0, model, schedule, name, 'NonLinearSolver', nls);

% We focus on the first 100 steps
n = 100;
ind = 1:n;
scheduleRef = problem.SimulatorSetup.schedule;
scheduleRef.step.val     = schedule.step.val(ind);
scheduleRef.step.control = schedule.step.control(ind);
Wref     = scheduleRef.control.W;
dt       = scheduleRef.step.val;
nstep    = numel(scheduleRef.step.val);

scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(n));

% For the training run, we create a schedule with varying controls
perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
rng(0)
scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .4, 'perturbStep', perturbStep);
% Pack the simulation problem with the initial state, model, and schedule
for i = 1:length(scheduleRef.control) 
    scheduleRef.control(i).bc = schedule.control(1).bc;
end

%%
problem.SimulatorSetup.schedule = scheduleRef;
simulatePackedProblem(problem, 'restartStep', 1);

[wsRef, statesRef] = getPackedSimulatorOutput(problem);
modelRef    = problem.SimulatorSetup.model;
%% Making Coarse model
x = 11;
y = 11;

map = modelRef.G.cells.indexMap;

% Nummerer indexMap unikt
modelRef.G.cells.indexMap = (1:numel(modelRef.G.cells.indexMap))';

blockIx = partitionUI(modelRef.G, [x, y]);
blockIx = processPartition(modelRef.G, blockIx);
blockIx = compressPartition(blockIx);
% Perform a simple upscaling to obtain a coarse model
modelRef.G.cells.indexMap = map;
modelCoarse = upscaleModelTPFA(modelRef, blockIx);
modelCoarse.AutoDiffBackend = AutoDiffBackend();
% We want to include rel-perm scaling as tunabale parameters, so include
% these for the coarse model. These parameters have no effect for the
% initial coarse model (they are set equal to the ones given by the
% rel-perm curves).

% Antall rader du ønsker
numRows = 3;

pts = modelCoarse.fluid.krPts;

scaling = {'SGL',   pts.g(1), 'SGCR', pts.g(2), 'SGU', pts.g(3), ...
            'SOGCR', pts.og(2), 'KRG',  pts.g(4), 'KRO', pts.og(4)};
modelCoarse = imposeRelpermScaling(modelCoarse, scaling{:});
modelCoarse.toleranceCNV = 1e-6;  % tighter tolerance to improve gradient accuracy

% % Utvid g og og til flere rader
modelCoarse.fluid.krPts.g = repmat(pts.g, numRows, 1);
modelCoarse.fluid.krPts.og = repmat(pts.og, numRows, 1);


%% Simulate initial upscaled coarse model and compare to reference
stateCoarse0   = upscaleState(modelCoarse, modelRef, state0);
% Changing cstatus which is which of the well cells that are active or not
stateCoarse0.wellSol.cstatus = logical([1, zeros(1, 7)]);
scheduleCoarse = upscaleSchedule(modelCoarse, scheduleRef);

% Changing cstatus which is which of the well cells that are active or not
% for i = 1:length(scheduleCoarse.control)
%     scheduleCoarse.control(i).W.cstatus = [true;false];
% end

%Validating the schedule:
scheduleCoarse = validateSchedule(modelCoarse, scheduleCoarse);

% Manually setting upper left to obtain symmetry
modelCoarse.rock.regions.saturation((6 * 11 + 2)) = 3;

%Manually moving the well one down...
% for i = 1:length(scheduleCoarse.control)
%     scheduleCoarse.control(i).W.cells = 61;
% end

% Manually setting the well to only be the bottom one:
%scheduleCoarse.control(1).W = addWell([], modelCoarse.G, modelCoarse.rock, ...
%    scheduleCoarse.control(1).W.cells(1), 'Type', 'rate', 'Val', scheduleCoarse.control(1).W.val, 'Radius', 0.1, 'name', 'W1');

%scheduleCoarse.control(1).W.compi = [1, 0];

%% Plot reference and coarse model grids with wells
plotWell2d(modelRef, modelCoarse, scheduleRef, scheduleCoarse)

%% Simulate
% Validating the model:
modelCoarse = modelCoarse.validateModel();
[wsCoarse, statesCoarse] = simulateScheduleAD(stateCoarse0, modelCoarse, scheduleCoarse);

plotWellSols({wsRef, wsCoarse}, ...
    {scheduleRef.step.val, scheduleCoarse.step.val},...
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

weighting  = {'OilRateWeight',  1/(sum(abs(qOs_ref))), ...         %Gas (earlier water): hydrogen       #Change so that this is also dependent on number of time steps (different cycles)
              'WaterRateWeight',    1/(sum(abs(qGs_ref))),...     %%Oil: water 
              'BHPWeight',        1/(max(abs(bhp_ref)))};


% make handle          
mismatchFn = @(model, states, schedule, states_ref, compDer, tstep, state) ...
    matchObservedH2W_changeBC(model, states, schedule, states_ref,...
                   'computePartials', compDer, 'tstep', tstep, weighting{:},...
                   'state', state, 'from_states', false);

%% Model calibration Quasi-Newton
pvec = getScaledParameterVector(setup_init, parameters);
objh = @(p) evaluateMatch(p, mismatchFn, setup_init, parameters, statesRef);
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

%% Compare convergence history:

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
plotWellSols({wsRef, wsCoarse, wellSols_opt1, wellSols_opt2}, ...
              repmat({scheduleRef.step.val}, 1, 4), ...
              'datasetnames', {'reference','coarse initial','coarse tuned, BFGS', 'coarse tuned (L-M)'});

%% Compare reference, initial coarse and optimizes coarse for a different schedule 
rng(100);
W = schedule.control.W;

dt = schedule.step.val(100:200);

s_new = perturbedSimpleSchedule(dt, 'W', Wref, ...
    'pressureFac', .01, 'rateFac', .2, 'perturbStep', ones(numel(dt),1));
ws_new1 = simulateScheduleAD(problem.SimulatorSetup.state0, modelRef, s_new, ...
    'NonLinearSolver', problem.SimulatorSetup.NonLinearSolver);

w_opt1 = setup_opt1.schedule.control(1).W;
w_opt2 = setup_opt2.schedule.control(1).W;
w_old = scheduleCoarse.control(1).W;   
s_opt1 = simpleSchedule(dt, 'W', w_opt1);
s_opt2 = simpleSchedule(dt, 'W', w_opt2);
s_old = simpleSchedule(dt, 'W', w_old);
%Setting controls:
for kw = 1:numel(Wref)
    s_opt1.control.W(kw).val = s_new.control.W(kw).val;
    s_opt2.control.W(kw).val = s_new.control.W(kw).val;
    s_old.control.W(kw).val = s_new.control.W(kw).val;
end

ws_coarse_new1 = simulateScheduleAD(setup_opt1.state0, setup_opt1.model, s_opt1);
ws_coarse_new2 = simulateScheduleAD(setup_opt2.state0, setup_opt2.model, s_opt2);
ws_coarse_old = simulateScheduleAD(stateCoarse0, modelCoarse, s_old);

%%

plotWellSols({ws_new1,  ws_coarse_old, ws_coarse_new1, ws_coarse_new2}, ...
    repmat({dt}, 1, 4), ...
    'datasetnames', {'reference','coarse initial','coarse tuned BFGS', 'coarse tuned LM'});

%% Plot solutions isolated
plotSolutions(ws_coarse_old, ws_coarse_new2)
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

function [description, options, state0, model, schedule, plotOptions] = H2_illustration_storage_example_cartesian_grid(deck,varargin)
% This example simulates the injection and behavior of hydrogen (H₂) in a 2D saline aquifer 
% using the black-oil model. The aquifer has a dome-shaped structure defined by the function 
% F(x) = σ + r*sin(π*x), with parameters σ = 25 and r = 5. 
% 
% The domain is 50 m × 50 m, with caprock and bedrock layers providing containment through 
% high entry (capillary) pressure and low permeability. 
% 
% A single well is placed at the top of the trap, in the uppermost cell beneath the caprock. 
% We apply fixed hydrostatic pressure at the lateral boundaries and enforce no-flux conditions 
% at the top and bottom boundaries.
%
% Relative permeability and capillary pressure follow the Brooks-Corey model, with varying 
% residual saturations for different rock types.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}
% Step 1: Test Case Description and Options
%---------------------------------------------------------------------%
% Description of the conceptual model for hydrogen storage with multiple
% injection/production (Inj/Prod) cycles.

description = 'Conceptual model for Hydrogen storage with multiple Inj/Prod cycles';

% Time units
K0    = 273.15 * Kelvin;  % Absolute temperature offset
% Optional input arguments
options = struct( ...
    'rateCharge'   , 18 * kilogram/day   , ... % Hydrogen injection rate during charging
    'rateIdle'     , 0.0 * kilogram/day  , ... % No injection during idle periods
    'rateCushion'  , 10 * kilogram/day   , ... % Cushion gas injection rate (H₂)
    'rateDischarge', 18 * kilogram/day   , ... % Hydrogen production rate during discharge
    'bhp'          , 35.0 * barsa        , ... % Bottom hole pressure during production
    'tempCharge'   , K0 + 40 * Kelvin    , ... % Injection temperature during charging
    'tempDischarge', K0 + 40 * Kelvin    , ... % Production temperature during discharge
    'tempCushion'  , K0 + 40 * Kelvin    , ... % Temperature for cushion gas phase
    'timeCushion'  , 90 * day            , ... % Duration of cushion gas phase
    'timeCharge'   , 30 * day            , ... % Duration of hydrogen charging
    'timeIdle'     , 10 * day            , ... % Idle period between charge/discharge cycles
    'timeShut'     , 30 * day            , ... % Shut-in period (no activity)
    'timeDischarge', 30 * day            , ... % Duration of hydrogen production (discharge)
    'dtCharge'     , 8.4 * hour          , ... % Timestep during charging
    'dtCushion'    , 8.4 * hour          , ... % Timestep during cushion phase
    'dtIdle'       , 8.4 * hour          , ... % Timestep during idle phase
    'dtShut'       , 8.4 * hour          , ... % Timestep during shut-in phase
    'dtDischarge'  , 8.4 * hour          , ... % Timestep during discharge
    'numCycles'    , 10                  , ... % Number of injection/production cycles
    'chargeOnly'   , 0                   , ... % Simulate only charging period
    'cushionOnly'  , 0                   , ... % Simulate only cushion gas phase
    'dischargeOnly', 0                   , ... % Simulate only discharge period
    'useGroupCtrl' , false               , ... % Group control for wells (optional)
    'initPres'     , 37 * barsa          , ... % Initial reservoir pressure
    'initTemp'     , K0 + 40             , ... % Initial reservoir temperature
    'initSat'      ,[1 0]                , ... % Initial reservoir saturation
    'use_bc'       , true                , ... % Use boundary conditions
    'use_cushion'  , true                , ... % Include cushion gas in the simulation
    'use_bhp'      , false                 ... % Use bottom hole pressure control
);

% Process optional input arguments
[options, fullSetup, ~] = processTestCaseInput(mfilename, ...
options, description, varargin{:});  % Process test case inputs
options = checkOptions(options);          % Check the validity of options

if ~fullSetup
    return;  % If setup is incomplete, exit early
end
% Merge optional input arguments with the existing options
options = merge_options(options, varargin{:});
% If the number of output arguments is less than or equal to 2, return early
if nargout <= 2
    return;
end

% Define module dependencies for the simulation
require ad-core ad-props ad-blackoil spe10 upr

% Grid setup
% Generate a constraint path along the x-axis and define the grid shape
x = linspace(0.0, 1.0);       % Generate linearly spaced points for x-axis
y = 25 + 5 * sin(pi * x);     % Define y-axis using a sinusoidal function for dome shape
w = {[50 .* x', y']};         % Combine x and y to form 2D points for grid constraints

% Grid size scaling and dimensions
gS = [0.5, 0.5];              % Scaling factors for grid
pdims = [50, 50];             % Grid dimensions in x and y directions

% Generate two versions of a composite PEBI grid
% G1 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
%     'interpolateCC', true, 'protLayer', false);  % Without protection layer

% G2 = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, ...
%     'interpolateCC', true, 'protLayer', true);   % With protection layer

% Compute geometry for the grid
G1 = cartGrid([100,100],[50,50]);
G = computeGeometry(G1);

  
% Slight adjustment to the y-axis values for accuracy
y = y + 0.1;  % Adjust y-values to ensure correct positioning above the grid
% Create a polygon representation of the line
linePolygon = [[0, 50]; [50 .* x', y']; [50, 50]];  % Define the polygon that represents the dome structure

% Define the cell vertices from the grid
cellVertices = G.cells.centroids;  % Extract centroids of the grid cells

% Use inpolygon to find which cells are above the line
aboveLineMask = inpolygon(cellVertices(:, 1), cellVertices(:, 2), linePolygon(:, 1), linePolygon(:, 2)); 

% Create rock properties
rock = makeRock(G, [10 * milli * darcy, 10 * milli * darcy], 0.25);  % Define rock permeability and porosity
caprock = aboveLineMask;  % Assign mask for caprock
bedrock = find(G.cells.centroids(:, 2) < 5);  % Identify bedrock cells based on their centroid positions

% Set permeability for caprock and bedrock
rock.perm(caprock, :) = 1.0e-04 * milli * darcy;  % Caprock permeability
rock.perm(bedrock, :) = 1.0e-02 * milli * darcy;  % Bedrock permeability

% Set porosity for caprock and bedrock
rock.poro(caprock) = 0.1;  % Caprock porosity
rock.poro(bedrock) = 0.1;  % Bedrock porosity

% Convert deck units for simulation
deck = convertDeckUnits(deck);
% Initialize activity and saturation numbers in the deck
deck.GRID.ACTNUM = ones(G.cells.num, 1);  % Set active cells
deck.REGIONS.SATNUM = ones(G.cells.num, 1);  % Set saturation numbers to 1

% Update saturation numbers for specific regions
deck.REGIONS.SATNUM(bedrock) = 2;  % Bedrock saturation number
deck.REGIONS.SATNUM(caprock) = 3;  % Caprock saturation number
rock.regions.saturation = deck.REGIONS.SATNUM;  % Assign saturation regions to rock properties
% We update indexmap and rock in G
G.cells.indexMap = rock.regions.saturation;
deck.GRID.PORO=rock.poro;
deck.GRID.PERMX=rock.perm(:,1);
deck.GRID.PERMY=rock.perm(:,2);
deck.GRID = rmfield(deck.GRID,'PERMZ');
% Initialize Eclipse problem for AD (Automatic Differentiation)
[~, model, ~] = initEclipseProblemAD(deck,'getSchedule',false,'G',G,'getInitialState', false);  

% Extract fluid and input data from the model
fluid = model.fluid;  
deck = model.inputdata;  


% Reset gravity for the model
gravity reset on;

% Create the black-oil model with gravity effects and input data
model = GenericBlackOilModel(G, rock, fluid, ...
    'water', false, 'disgas', true, 'gravity', [0 -norm(gravity)], 'inputdata', deck);

% Set up the simulation schedule based on the grid, rock properties, fluid model, and options
schedule = setUpSchedule(G, rock, fluid, options);

% Set up the initial state for the simulation using the first control well from the schedule
state0 = setUpInitialState(model, schedule.control(1).W, options);

% Define plotting options for visualization
plotOptions = {'View'              , [0,0]         , ...
    'PlotBoxAspectRatio', [1,1,0.25]    , ...
    'Projection'        , 'orthographic', ...
    'Size'              , [800, 300]    };
% Optional: Set title for the simulation run in the deck (commented out)
% deck.RUNSPEC.TITLE = 'H2_illustration_storage';

% Optional: Convert the model to a deck format (commented out)
% deck_new = model2Deck(model, schedule, 'deck', deck);

end

function W = setUpWells(G, rock, fluid, options)
    % setUpWells - Sets up the wells for the simulation
    % 
    % Syntax: W = setUpWells(G, rock, fluid, options)
    %
    % Inputs:
    %   G      - Grid structure
    %   rock   - Rock properties structure
    %   fluid  - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   W      - Well structure containing the defined wells

    % Initialize the well array
    W = [];  
    
    % Identify cells for the production well based on specified coordinates
    wc = find(abs(G.cells.centroids(:, 1) - 25.25) < 0.1 & ...
              abs(G.cells.centroids(:, 2)) > 25.5 & ...
              abs(G.cells.centroids(:, 2)) < 29.5);
    
    % Add a production well at the identified cells with specified properties
    W = addWell(W, G, rock, wc, ...
                'Name', 'Prod', ...                       % Well name
                'Radius', 5 * centi * meter, ...        % Well radius
                'Type', 'rate', ...                      % Well type (rate control)
                'Val', options.rateCharge, ...           % Production rate
                'Compi', [0, 1]);                        % Component indices

    % Set well groups if group control is enabled
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'});         % Assign groups for injection and production
    end
end


function bc = setUpBc(G, rock, fluid, options)
    % setUpBc - Sets up boundary conditions for the simulation.
    %
    % Syntax: bc = setUpBc(G, rock, fluid, options)
    %
    % Inputs:
    %   G      - Grid structure
    %   rock   - Rock properties structure
    %   fluid  - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   bc     - Boundary condition structure

    % Initialize saturation and fluid density
    sat = options.initSat;
    omega = fluid.rhoOS;  % Water saturation density

    % Determine gravity vector based on grid dimensions
    if G.griddim < 3
        grav_ = [0 -norm(gravity)];  % 2D case
    else
        grav_ = gravity;  % 3D case
    end

    % Identify boundary faces
    f = boundaryFaces(G);
    f1 = f(abs(G.faces.centroids(f, 1)) < eps);     % Left boundary
    f2 = f(abs(G.faces.centroids(f, 1) - 50) < eps); % Right boundary

    % Calculate pressure difference for left boundary
    dx1 = bsxfun(@minus, G.faces.centroids(f1, :), 0);
    dp1 = omega .* (dx1 * reshape(grav_, [], 1));  % Pressure change due to gravity
    p0 = options.initPres;                           % Initial pressure
    pcmax = p0 + 3 * barsa();                       % Maximum pressure condition
    pressure1 = pcmax + dp1;                           % Total pressure for left boundary

    % Add boundary condition for left boundary
    bc = addBC([], f1, 'pressure', pressure1, 'sat', sat);

    % Calculate pressure difference for right boundary
    dx2 = bsxfun(@minus, G.faces.centroids(f2, :), 0);
    dp2 = omega .* (dx2 * reshape(grav_, [], 1));  % Pressure change due to gravity
    pressure2 = pcmax + dp2;                       % Total pressure for right boundary

    % Add boundary condition for right boundary
    bc = addBC(bc, f2, 'pressure', pressure2, 'sat', sat);
end



% function schedule = setUpSchedule(G0, rock, fluid, options)
%         
%     W = setUpWells(G0, rock, fluid, options);     
%     if options.use_cushion
%        W(1).type     = 'rate';
%        W(1).name     = 'cushion';
%        W(1).val      = options.rateCushion;
%        W(1).T        = options.tempCushion;
%        W(1).sign     = 1;
%     
%        dtCushions       = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
%        rateCushion = options.rateCushion.*dtCushions./max(dtCushions);
%        for i = 1:9    
%            dtCushion      = dtCushions(i);
%            W(1).val      = rateCushion(10);
%            scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
%        end
% 
%        dtCushion       = dtCushions(10:end-10);
%        W(1).val      = rateCushion(10);
%        scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);
% 
%        
%        for i = 1:9    
%            dtCushion      = dtCushions(end-9+i);
%            W(1).val      = rateCushion(end-9+i);
%            scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
%        end
%     end
%     W = setUpWells(G0, rock, fluid, options);    
%     W(1).type     = 'rate';
%     W(1).name     = 'charge';    
%     W(1).val      = options.rateCharge;
%     W(1).T        = options.tempCharge;
%     W(1).sign     = 1;
%     
%     dtCharges       = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
%     rateCharge = options.rateCharge.*dtCharges./max(dtCharges);
%     for i = 1:9    
%         dtCharge      = dtCharges(i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
%     end
% 
%     dtCharge       = dtCharges(10:end-10);
%     W(1).val      = rateCharge(10);
%     scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);
% 
%        
%     for i = 1:9    
%         dtCharge      = dtCharges(end-9+i);
%         W(1).val      = rateCharge(10);
%         scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
%     end
%     
%     
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
% 
%     W(1).type     = 'rate';
%     W(1).val      = options.rateIdle;
%     W(1).name     = 'shut';        
%     W(1).T        = options.tempCushion;
%     W(1).sign     = -1;
%     
%     dtIdle       = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
%     scheduleIdle = simpleSchedule(dtIdle, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleIdle.groups = groups;
%     end
% 
%     dtShut       = rampupTimestepsEnds(options.timeShut, options.dtShut);
%     scheduleShut = simpleSchedule(dtShut, 'W', W);
% 
%     W(1).type     = 'rate';
%     W(1).name     = 'discharge';    
%     W(1).val      = -options.rateDischarge;
%     W(1).sign     = -1;
%     
%     dtDischarge       = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
%     scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);
%     if options.useGroupCtrl
%         groups = [];
%         scheduleCharge.groups = groups;
%     end
%     
%     if options.chargeOnly
%         schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
%     elseif options.dischargeOnly
%         schedule = scheduleDischarge;
%     elseif options.cushionOnly
%         schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
%     else 
%         if options.use_cushion
%            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%            schedule = repmat({schedule}, 1, options.numCycles);
%            schedule = combineSchedules(scheduleCushions{:}, scheduleShut, schedule{:}, 'makeConsistent', false);
%         else
%             schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
%             schedule = repmat({schedule}, 1, options.numCycles);
%             schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
%         end
%     end
% 
%         
%     if options.use_bc    
%         bc = setUpBc(G0,rock,fluid,options);        
%         for i = 1:numel(schedule.control)        
%             schedule.control(i).bc = bc;
%         end
%     end
% 
% end
function schedule = setUpSchedule(G0, rock, fluid, options)
    % setUpSchedule - Sets up the well schedule for the reservoir simulation.
    %
    % Syntax: schedule = setUpSchedule(G0, rock, fluid, options)
    %
    % Inputs:
    %   G0      - Grid structure
    %   rock    - Rock properties structure
    %   fluid   - Fluid properties structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   schedule - Schedule structure containing well control information

    % Set up wells
    W = setUpWells(G0, rock, fluid, options); 

    % Define cushion schedule if applicable
    if options.use_cushion
        W(1).type = 'rate';
        W(1).name = 'cushion';
        W(1).val = options.rateCushion;
        W(1).T = options.tempCushion;
        W(1).sign = 1;

        % Ramp up cushion rates
        dtCushions = rampupTimestepsEnds(options.timeCushion, options.dtCushion);
        rateCushion = options.rateCushion .* dtCushions ./ max(dtCushions);
        
        % Create cushion schedules
        for i = 1:9    
            dtCushion = dtCushions(i);
            W(1).val = rateCushion(10);  % Use the tenth value for the first nine steps
            scheduleCushions{i} = simpleSchedule(dtCushion, 'W', W);
        end

        % Handle remaining cushion schedules
        dtCushion = dtCushions(10:end-10);
        W(1).val = rateCushion(10);
        scheduleCushions{10} = simpleSchedule(dtCushion, 'W', W);

        for i = 1:9    
            dtCushion = dtCushions(end-9+i);
            W(1).val = rateCushion(end-9+i);
            scheduleCushions{i+10} = simpleSchedule(dtCushion, 'W', W);
        end
    end

    % Define charge schedule
    W = setUpWells(G0, rock, fluid, options); 
    W(1).type = 'rate';
    W(1).name = 'charge';    
    W(1).val = options.rateCharge;
    W(1).T = options.tempCharge;
    W(1).sign = 1;

    % Ramp up charge rates
    dtCharges = rampupTimestepsEnds(options.timeCharge, options.dtCharge);
    rateCharge = options.rateCharge .* dtCharges ./ max(dtCharges);
    
    % Create charge schedules
    for i = 1:9    
        dtCharge = dtCharges(i);
        W(1).val = rateCharge(10);  % Use the tenth value for the first nine steps
        scheduleCharges{i} = simpleSchedule(dtCharge, 'W', W);
    end

    % Handle remaining charge schedules
    dtCharge = dtCharges(10:end-10);
    W(1).val = rateCharge(10);
    scheduleCharges{10} = simpleSchedule(dtCharge, 'W', W);

    for i = 1:9    
        dtCharge = dtCharges(end-9+i);
        W(1).val = rateCharge(end-9+i);
        scheduleCharges{i+10} = simpleSchedule(dtCharge, 'W', W);
    end

    % Set up idle schedule
    W(1).type = 'rate';
    W(1).val = options.rateIdle;
    W(1).name = 'shut';        
    W(1).T = options.tempCushion;
    W(1).sign = -1;

    dtIdle = rampupTimestepsEnds(options.timeIdle, options.dtIdle);
    scheduleIdle = simpleSchedule(dtIdle, 'W', W);

    % Set up shut schedule
    dtShut = rampupTimestepsEnds(options.timeShut, options.dtShut);
    scheduleShut = simpleSchedule(dtShut, 'W', W);

    % Define discharge schedule
    W(1).type = 'rate';
    W(1).name = 'discharge';    
    W(1).val = -options.rateDischarge;
    W(1).sign = -1;

    dtDischarge = rampupTimestepsEnds(options.timeDischarge, options.dtDischarge);
    scheduleDischarge = simpleSchedule(dtDischarge, 'W', W);

    % Combine schedules based on options
    if options.chargeOnly
        schedule = combineSchedules(scheduleCharges{:}, 'makeConsistent', false);
    elseif options.dischargeOnly
        schedule = scheduleDischarge;
    elseif options.cushionOnly
        schedule = combineSchedules(scheduleCushions{:}, 'makeConsistent', false);
    else 
        if options.use_cushion
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(scheduleCushions{:}, scheduleShut, schedule{:}, 'makeConsistent', false);
        else
            schedule = combineSchedules(scheduleCharges{:}, scheduleIdle, scheduleDischarge, 'makeConsistent', false);
            schedule = repmat({schedule}, 1, options.numCycles);
            schedule = combineSchedules(schedule{:}, 'makeConsistent', false);  
        end
    end

    % Set up boundary conditions if specified
    if options.use_bc    
        bc = setUpBc(G0, rock, fluid, options);        
        for i = 1:numel(schedule.control)        
            schedule.control(i).bc = bc;  % Assign boundary conditions to each control
        end
    end
end

function state0 = setUpInitialState(model, W, options)
    % setUpInitialState - Initializes the state of the reservoir simulation.
    %
    % Syntax: state0 = setUpInitialState(model, W, options)
    %
    % Inputs:
    %   model  - Simulation model structure containing grid and rock properties
    %   W      - Well structure
    %   options - Simulation options structure
    %
    % Outputs:
    %   state0 - Initial state structure containing pressure and saturation

    % Initialize reservoir solution with specified initial pressure and default saturations
    state0 = initResSol(model.G, options.initPres, options.initSat);  % [1, 0] represents initial saturations

    % Initialize residual saturations (set to zero initially)
    state0.rs = 0 .* state0.pressure;

    % Uncomment and modify the following lines to adjust saturations based on bedrock properties
    % bedrock = find(model.G.cells.centroids(:, 2) < 15);
    % state0.s(bedrock, 1) = 1 - (model.rock.poro(bedrock, 1) .* model.G.cells.centroids(bedrock, 2) ./ 5);
    % state0.s(bedrock, 2) = (model.rock.poro(bedrock, 1) .* model.G.cells.centroids(bedrock, 2) ./ 5);

    % Initialize well solutions
    wellSol = initWellSolAD(W, model, state0);
    wellSol.bhp = options.initPres;  % Set initial bottom hole pressure for wells
    state0.wellSol = wellSol;         % Assign well solutions to the state

end


%-------------------------------------------------------------------------%
function options = checkOptions(options)
    
    assert(~(options.chargeOnly && options.dischargeOnly), ...
        'Cannot simulate only charge and only discharge at the same time');
    
end

function dT = rampupTimestepsEnds(time, dt, n)
% Create timesteps that ramp up geometrically
%
% SYNOPSIS:
%   dT = rampupTimesteps(1*year, 30*day)
%   dT = rampupTimesteps(1*year, 30*day, 5)
%
% DESCRIPTION:
%   This function generates a timestep sequence for a given total time
%   interval that increases geometrically until it reaches some target
%   timestep. The rest of the interval is then divided into a number of
%   target timesteps.
%
% REQUIRED PARAMETERS:
%   time   - The total simulation time so that sum(dt) = time
%
%   dt     - Target timestep after initial ramp-up
%
%   n      - (OPTIONAL) Number of rampup steps. Defaults to 8.
%
% RETURNS:
%   dt     - Array of timesteps.
%
% NOTE:
%   The final timestep may be shorter than dt in order to exactly reach T.
%

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if nargin < 3
        n = 8;
    end
    if time == 0
        dT = [];
        return
    end
    % Initial geometric series
    dt_init = (dt./2.^[n n:-1:1])';
    cs_time = cumsum(dt_init);
    if any(cs_time > time)
        dt_init = dt_init(cs_time < time);
    end
    
    % Remaining time that must be discretized
    dt_left = time - sum(2.*dt_init);
    % Even steps
    dt_rem = repmat(dt, floor(dt_left/dt), 1);
    % Final ministep if present
    dt_final = time - sum(2.*dt_init) - sum(dt_rem);
    % Less than to account for rounding errors leading to a very small
    % negative time-step.
    if dt_final <= 0
        dt_final = [];
    end
       
    if dt_final >= dt_init(1)
        dt_final = dt_init(1);
    end
    % Combined timesteps
    dT = [dt_init; dt_rem;sort(dt_init,'descend'); dt_final];
end

function obj = matchObservedH2W_changeBC(model, states, schedule, observed, varargin)
% Compute mismatch-function 
%previously named 'matchObservedOW' for oil water mismatch-function. Now we
%have gas water 

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
opt     = struct('WaterRateWeight',     [] , ...    %WaterRateWeight er min WaterRateWeight (notation O)
                 'OilRateWeight',       [] , ...    %OilRateWeight er min HydrogengasRateWeight (notation G)
                 'BHPWeight',           [] , ...
                 'ComputePartials',     false, ...
                 'tStep' ,              [], ...
                 'state',               [],...
                 'from_states',         true,...% can be false for generic models
                 'matchOnlyProducers',  false, ...
                 'mismatchSum',         true, ...
                 'accumulateWells',       [], ...
                 'accumulateTypes',       []);
             
opt     = merge_options(opt, varargin{:});

dts   = schedule.step.val;
totTime = sum(dts);

tSteps = opt.tStep;
if isempty(tSteps) %do all
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    numSteps = 1;
    dts = dts(opt.tStep);
end


obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    sol_obs = observed{tSteps(step)};
    nw      = numel(sol_obs.wellSol);
    if opt.matchOnlyProducers
        matchCases = (vertcat(sol.sign) < 0);
    else
        matchCases = true(nw,1);
    end
    qOs_obs = vertcatIfPresent(sol_obs.wellSol, 'qOs', nw);
    qGs_obs = vertcatIfPresent(sol_obs.wellSol, 'qGs', nw);
    bhp_obs = vertcatIfPresent(sol_obs.wellSol, 'bhp', nw);
  
    status_obs = vertcat(sol_obs.wellSol.status);
    
    [ww, wg, wp] = getWeights(qOs_obs, qGs_obs, bhp_obs, opt);
    
    if opt.ComputePartials
        if(opt.from_states)
            init=true;
            state = model.getStateAD( states{tSteps(step)}, init);
        else
            state = opt.state;
        end
        qGs = model.FacilityModel.getProp(state,'qGs');
        qOs = model.FacilityModel.getProp(state,'qOs');
        bhp = model.FacilityModel.getProp(state,'bhp');
        
        assert(not(isnumeric(qOs))); 
        status = vertcat(state.wellSol.status);
     else
        state = states{tSteps(step)};
        [qOs, qGs, bhp] = deal( vertcatIfPresent(state.wellSol, 'qOs', nw), ...
                                vertcatIfPresent(state.wellSol, 'qGs', nw), ...
                                vertcatIfPresent(state.wellSol, 'bhp', nw) );
       assert(isnumeric(qOs));
       status = vertcat(state.wellSol.status);
    end
     
    if ~all(status) || ~all(status_obs) 
        [bhp, bhp_obs] = expandToFull(bhp, bhp_obs, status, status_obs, true);
        [qOs, qOs_obs] = expandToFull(qOs, qOs_obs, status, status_obs, false);
        [qGs, qGs_obs] = expandToFull(qGs, qGs_obs, status, status_obs, false);
    end
    dt = dts(step);
    if opt.mismatchSum %building functional
        obj{step} = (dt/(totTime*nnz(matchCases)))*sum( ...
                        (ww.*matchCases.*(qOs-qOs_obs)).^2 + ...     %*ww
                        (wg.*matchCases.*(qGs-qGs_obs)).^2 + ...     %*wo  
                        (wp.*matchCases.*(bhp-bhp_obs)).^2);% + ...        %*wp
                        %(1.*matchCases.*(flux - flux_obs)./flux_obs).^2);   %changed here: flux along boundaries 
    else
        % output summands f_i^2 
        fac = dt/(totTime*nnz(matchCases));
        mm  = {fac*(ww.*matchCases.*((qOs-qOs_obs))).^2, ...           %*ww
               fac*(wg.*matchCases.*((qGs-qGs_obs))).^2, ...           %*wo
               fac*(wp.*matchCases.*((bhp-bhp_obs))).^2}; %, ...              %*wp
                %fac*(matchCases.*((flux - flux_obs)./flux_obs)).^2};    %changed here: fluix along boundaries
        
        if isempty(opt.accumulateTypes)
            tmp = mm;
        else
            % sum squares of qGs/qOs/bhp
            pt = opt.accumulateTypes;
            tmp = num2cell(zeros(1, max(pt)));
            for k = 1:3
                if pt(k)>0
                    tmp{pt(k)} = tmp{pt(k)} + mm{k};
                end
            end
        end
        if ~isempty(opt.accumulateWells)
            % sum squares of values for wells (use sparse mult)
            pw  = opt.accumulateWells;
            M   = sparse(pw(pw>0), find(pw), 1);
            tmp = applyFunction(@(x)M*x, tmp);
        end
        obj{step} = vertcat(tmp{:});
    end
end

end
%--------------------------------------------------------------------------
function v = vertcatIfPresent(sol, fn, nw)
if isfield(sol, fn)
    v = vertcat(sol.(fn));
    assert(numel(v)==nw);
    v = v(vertcat(sol.status));
else
    v = zeros(nnz(sol.status),1);
end
end

%--------------------------------------------------------------------------

function [v, v_obs] = expandToFull(v, v_obs, status, status_obs, setToZero)
tmp = zeros(size(status));
if isa(v, 'ADI')
    tmp = double2ADI(tmp, v);
end
tmp(status) = v;
v = tmp;
%
tmp = zeros(size(status));
tmp(status_obs) = v_obs;
v_obs = tmp;
if setToZero
    ix = status ~= status_obs;
    v(ix)     = 0;
    v_obs(ix) = 0;
end

end
%--------------------------------------------------------------------------

function  [ww, wo, wp] = getWeights(qWs, qOs, bhp, opt)
ww = opt.WaterRateWeight;
wo = opt.OilRateWeight;
wp = opt.BHPWeight;

rw = sum(abs(qWs)+abs(qOs));

if isempty(ww)
    % set to zero if all are zero
    if sum(abs(qWs))==0
        ww = 0;
    else
        ww = 1/rw;
    end
end

if isempty(wo)
    % set to zero if all are zero
    if sum(abs(qOs))==0
        wo = 0;
    else
        wo = 1/rw;
    end
end

if isempty(wp)
    % set to zero all are same
    dp = max(bhp)-min(bhp);
    if dp == 0
        wp = 0;
    else
        wp = 1/dp;
    end
end
end