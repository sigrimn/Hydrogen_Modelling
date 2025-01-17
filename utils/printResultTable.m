function result_table = printResultTable(results, wsRef_full)

[v1, p_opt1, history1, v2, p_opt2, history2, modelCoarse, modelRefPert] = loadResults(results);

%% Evaluerer resultater
nit_BFGS = length(history1.val);
nit_LM = length(history2.val);

err_BFGS = v1;
err_LM = v2;

[bhp_error_init, h2o_rate_error_init, h2_rate_error_init] = calculate_well_errors(wsRef_full, results.ws_init);
[bhp_error, h2o_rate_error, h2_rate_error] = calculate_well_errors(wsRef_full, results.ws_opt);

%% Setter sammen tabell

% Sett sammen data i en tabell
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

% Skriv ut tabellen
disp('Resultater for simuleringen:');
disp(result_table);