function params = setupParameters(setup_init, config)
    % Set up model parameters from config, matching
    % the model dimensions fo the network model netmod
    
    if isempty(config)
        warning("config is empty. No parameters constructed.")
    end
    params = [];
    for k = 1:size(config,1)
    if config{k, 2} == 0, continue, end
    params = addParameter(params, setup_init, ...
        'name',    config{k,1}, 'scaling', config{k,3}, ...
        'boxLims', config{k,4}, 'relativeLimits',config{k,5});
end
end