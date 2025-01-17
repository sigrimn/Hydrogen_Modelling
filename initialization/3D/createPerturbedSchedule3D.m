function finalSchedule = createPerturbedSchedule3D(ref, nbhp, ninj, nprod, nshutt, pressureFac, rateFac)


%% Optional input arguments: options corresponds to scenario 2 of the paper
K0 = 273.15 * Kelvin;
options = struct( ...
    'rateCharge',       6.14062 * 10^5 * meter^3/day/2,   ... % Charge rate (m³/day)
    'rateIdle',         0.0 * kilogram/day/2,             ... % Idle rate (kg/day)
    'rateCushion',      6.14062 * 10^5 * meter^3/day/2,   ... % Cushion (H₂) injection rate (m³/day)
    'rateDischarge',    6.14062 * 10^5 * meter^3/day/2,   ... % Discharge rate (m³/day)
    'bhp',              90.0 * barsa,                     ... % Production bottom hole pressure (BHP) in bar
    'tempCharge',       K0 + 60 * Kelvin,                 ... % Charging temperature (K)
    'tempDischarge',    K0 + 60 * Kelvin,                 ... % Discharging temperature (K)
    'tempCushion',      K0 + 60 * Kelvin,                 ... % Cushion temperature (K)
    'timeCushion',      120 * day,                        ... % Duration of cushioning phase (days)
    'timeCharge',       30 * day,                         ... % Duration of charging phase (days)
    'timeIdle',         15 * day,                         ... % Duration of idle phase (days)
    'timeShut',         30 * day,                         ... % Duration of shut phase (days)
    'timeDischarge',    30 * day,                         ... % Duration of discharging phase (days)
    'dtCharge',         1.0 * day,                        ... % Time step during charging (days)
    'dtCushion',        1.0 * day,                        ... % Time step during cushioning (H₂) (days)
    'dtIdle',           1.0 * day,                        ... % Time step during idle (days)
    'dtShut',           1.0 * day,                        ... % Time step during shut phase (days)
    'dtDischarge',      1.0 * day,                        ... % Time step during discharging (days)
    'numCycles',        20,                               ... % Total number of cycles (charging and discharging)
    'numCyclesCushions',6,                                ... % Number of cycles for cushion gas (H₂)
    'chargeOnly',       0,                                ... % Flag to simulate only the charging phase (0: No, 1: Yes)
    'cushionOnly',      0,                                ... % Flag to simulate only the cushioning phase (0: No, 1: Yes)
    'dischargeOnly',    0,                                ... % Flag to simulate only the discharging phase (0: No, 1: Yes)
    'useGroupCtrl',     false,                            ... % Flag to use group control (true/false)
    'initPres',         81.6 * barsa,                     ... % Initial pressure (bar)
    'initSat',          [1, 0],                           ... % Initial saturation
    'initTemp',         273.15 + 60,                      ... % Initial temperature (°C)
    'use_bc',           true,                             ... % Flag to use boundary conditions (true/false)
    'use_cushion',      true,                             ... % Flag to use cushion phase (true/false)
    'use_bhp',          true                              ... % Flag to use bottom hole pressure (true/false)
    );


%Initializing well
%Wref = setUpWells(G0, rock, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% bottom hole pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nbhp > 0
        ind = 1:nbhp;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref_BHP     = scheduleRef.control.W;

        %Wref_BHP = Wref;
        Wref_BHP(1).type = 'bhp';
        Wref_BHP(1).name = 'cushion';
        Wref_BHP(1).val = options.bhp;
        Wref_BHP(1).T = options.tempCushion;
        Wref_BHP(1).sign = 1;
        rng(0)

        dt       = scheduleRef.step.val;
        nstep    = numel(scheduleRef.step.val);
        
        scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(nbhp));
        
        % For the training run, we create a schedule with varying controls
        perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
        rng(0)
        scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref_BHP, ...
            'pressureFac', pressureFac, 'rateFac', rateFac, 'perturbStep', perturbStep);  %0.01 og 0.4 først!!!
        % Pack the simulation problem with the initial state, model, and schedule
        for i = 1:length(scheduleRef.control) 
            scheduleRef.control(i).bc = ref.control(1).bc;
        end

        schedule_bhp = scheduleRef;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% injection rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ninj > 0
        
        ind = 1:ninj;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref_INJ     = scheduleRef.control.W;

        Wref_INJ(1).type = 'rate';
        Wref_INJ(1).name = 'charge';
        Wref_INJ(1).val = options.rateCharge;
        Wref_INJ(1).T = options.tempCharge;
        Wref_INJ(1).sign = 1;
        rng(0)

        dt       = scheduleRef.step.val;
        nstep    = numel(scheduleRef.step.val);
        
        scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(ninj));
        
        % For the training run, we create a schedule with varying controls
        perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
        rng(0)
        scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref_INJ, ...           %endret her fra Wref_INJ !!!!
            'pressureFac', pressureFac, 'rateFac', rateFac, 'perturbStep', perturbStep);
        % Pack the simulation problem with the initial state, model, and schedule
        for i = 1:length(scheduleRef.control) 
            scheduleRef.control(i).bc = ref.control(1).bc;
        end
        schedule_inj = scheduleRef;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% production/discharge rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nprod > 0
        ind = 1:nprod;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref_PROD     = scheduleRef.control.W;

        Wref_PROD(1).type = 'rate';
        Wref_PROD(1).name = 'discharge';
        Wref_PROD(1).val = -options.rateDischarge;
        Wref_PROD(1).T = options.tempCharge;
        Wref_PROD(1).sign = -1;
        Wref_PROD(1).lims.bhp = 20*barsa;
        Wref_PROD(1).cstatus(2:end) = 0;
        rng(0)

        dt       = scheduleRef.step.val;
        nstep    = numel(scheduleRef.step.val);
        
        scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(nprod));
        
        % For the training run, we create a schedule with varying controls
        perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
        rng(0)
        scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref_PROD, ...
            'pressureFac', pressureFac, 'rateFac', rateFac, 'perturbStep', perturbStep);
        % Pack the simulation problem with the initial state, model, and schedule
        for i = 1:length(scheduleRef.control) 
            scheduleRef.control(i).bc = ref.control(1).bc;
        end
        schedule_prod= scheduleRef;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% shutting rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nshutt > 0
        ind = 1:nshutt;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref_SHUTT     = scheduleRef.control.W;

        Wref_SHUTT = scheduleRef.control.W;
        Wref_SHUTT(1).type = 'rate';
        Wref_SHUTT(1).name = 'charge';
        Wref_SHUTT(1).val = options.rateCharge;
        Wref_SHUTT(1).T = options.tempCharge;
        Wref_SHUTT(1).sign = 1;
        rng(0)



        dt       = scheduleRef.step.val;
        nstep    = numel(scheduleRef.step.val);
        
        scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(nshutt));
        
        % For the training run, we create a schedule with varying controls
        perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
        rng(0)
        scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref_SHUTT, ...
            'pressureFac', pressureFac, 'rateFac', rateFac, 'perturbStep', perturbStep);
        % Pack the simulation problem with the initial state, model, and schedule
        for i = 1:length(scheduleRef.control) 
            scheduleRef.control(i).bc = ref.control(1).bc;
        end

        schedule_bhp = scheduleRef;
        
    end

    %% COMBINE SCHEDULES

   % Combine schedules in the specified order
    %% COMBINE SCHEDULES

    % Initialize an empty cell array to store valid schedules
    combinedSchedules = {};

    % Add each schedule only if it exists (length > 0)
    if nbhp > 0 
        combinedSchedules = {combinedSchedules{:}, schedule_bhp}; 
    end

    if ninj > 0
        combinedSchedules = {combinedSchedules{:}, schedule_inj}; 
    end

    if nprod > 0
        combinedSchedules = {combinedSchedules{:}, schedule_prod}; 
    end

    if nshutt > 0
        combinedSchedules = {combinedSchedules{:}, schedule_shutt}; 
    end


    % Combine all valid schedules
    if ~isempty(combinedSchedules)
        % Use combineSchedules to merge all non-zero-length schedules
        finalSchedule = combineSchedules(combinedSchedules{:}, 'makeConsistent', false);
    else
        % If no valid schedules, return an empty schedule
        finalSchedule = struct();
        finalSchedule.step = struct('val', [], 'control', []);
        finalSchedule.control = struct();
    end
end 

function W = setUpWells(G, rock, options)

    %% We also reset the well in the highest point: we reset Well coordinates (wc) and radii (r)
    % ps this is slightly different from the benchmark
    wc = [5891; 9612; 13333; 17054; 20775; 24496; 28217; 31938; 35659; 39380; 43101];
    r = [0.9479; 0.7354; 0.7338; 2.2337; 2.2337; 2.6778; 2.6722; 5.1490; 5.1491; 0.7834; 0.7834];
    
    %% Add production wells with specified parameters
    W = addWell([], G, rock, wc, 'Name', 'Prod', ...
                'Radius', r, 'type', 'rate', ...
                'val', options.rateCharge, 'compi', [0, 1]);
    
    %% Set groups if group control is used
    if options.useGroupCtrl
        [W.group] = deal({'Inj', 'Prod'});
    end

end


