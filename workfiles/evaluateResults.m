%% Startup
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

%% Paths to work/util files
addpath(genpath('/Users/sigridmarianese/Documents/ProsjektOppgave'))
addpath(genpath('/Applications/MRST/modules/H2store'))

% Fetching 2d or 3d reference model for comparison
name3D = 'H2_STORAGE_RS_UHS_BENCHMARK_BC_SHORT_PERIDOS_LOWRATES_SALT';
deck = readEclipseDeck('/Users/sigridmarianese/Documents/MATLAB_workfolder/HydrogenOptimization/UHS/UHS_benchmark/UHSDATA/UHS_BENCHMARK_RS.DATA');

% Prepare simulation parameters and initialize
problem = initialize3D(deck, name3D);
[wsRef_full, statesRef] = getPackedSimulatorOutput(problem);
schedule = problem.SimulatorSetup.schedule;

%% Fetching results for comparison
resultspath = "/Users/sigridmarianese/Documents/ProsjektOppgave/3Dresults_rate_perturbation/";
resultsname = "UHS_05_001.mat";

result_ = strcat(resultspath, resultsname);
disp(resultspath + resultsname);
load(result_);

%% Evaluating results
nit_BFGS = length(results.history1.val);
nit_LM = length(results.history2.val);

err_BFGS = results.v1;
err_LM = results.v2;

[bhp_error_init, h2o_rate_error_init, h2_rate_error_init] = calculate_well_errors(wsRef_full, results.ws_init);
[bhp_error, h2o_rate_error, h2_rate_error] = calculate_well_errors(wsRef_full, results.ws_opt);

%% Assemble table

% Gathering data to table
result_table = table(...
    nit_BFGS, nit_LM, ...                                   % Iterasjoner
    err_BFGS, err_LM, ...                                   % Objektivfunksjonsfeil
    bhp_error_init, bhp_error, ...                          % BHP-feil (init og opt)
    h2o_rate_error_init, h2o_rate_error, ...                % H2O-rate feil (init og opt)
    h2_rate_error_init, h2_rate_error, ...                  % H2-rate feil (init og opt)
    'VariableNames', { ...
        'BFGS_Iter', 'LM_Iter', ...                         % Kolonnetitler for iterasjoner
        'Err_BFGS', 'Err_LM', ...                           % Kolonnetitler for feil
        'BHP_Error_Init', 'BHP_Error_Opt', ...              % Kolonnetitler for BHP-feil
        'H2O_Rate_Error_Init', 'H2O_Rate_Error_Opt', ...    % Kolonnetitler for H2O-rate feil
        'H2_Rate_Error_Init', 'H2_Rate_Error_Opt' ...       % Kolonnetitler for H2-rate feil
    });

% Display table
disp('Resultater for simuleringen:');
disp(result_table);

%% Save table
writetable(result_table, 'resultater.csv');



%%
% Plotting optimized results
plotWellSols({results.wsRefPert, results.wsCoarse, results.wellSols_opt1, results.wellSols_opt2}, ...
              repmat({results.results.scheduleRef.step.val}, 1, 4), ...
              'datasetnames', {'reference','coarse initial','coarse tuned (BFGS)', 'coarse tuned (LM)'});


%% Plotting full simulations
ind = 1:1400;
plotWellSols({results.wsRef_full(ind), results.wsCoarse(ind), results.ws_opt_BFGS(ind), results.ws_opt_LM(ind)}, ...
    repmat({schedule.step.val(ind)}, 1, 4), ...
    'datasetnames', {'reference', 'coarse init', 'coarse tuned (BFGS)', 'coarse tuned (LM)'});


%% Plotting history
figure
semilogy(abs(results.history1.val)', '-o', 'LineWidth', 2);
hold on
semilogy(abs(results.history2.val)', '--o', 'LineWidth', 2);
hold on
semilogy(abs(lm_static_history.val)', '-o', 'LineWidth', 2);
set(gca, 'Fontsize', 12); grid on
xlabel('Iteration'), ylabel('Residual')
legend({'LM', 'BFGS + LM', 'LM'})