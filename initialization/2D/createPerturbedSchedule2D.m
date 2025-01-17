function finalSchedule = createPerturbedSchedule2D(ref, nbhp, ninj, nprod, nshutt, pressureFac, rateFac)



% scheduleOrder: order of controls. 
% example: scheduleOrder = {'bhp', 'inj'} means the bhp control is first,
% inj control last


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% bottom hole pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nbhp > 0
        ind = 1:nbhp;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);

        Wref_BHP = scheduleRef.control.W;
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
        Wref     = scheduleRef.control.W;

        Wref_INJ = scheduleRef.control.W;
        Wref_INJ(1).type = 'rate';
        Wref_INJ(1).name = 'cushion';
        Wref_INJ(1).val = options.rateCushion;
        Wref_INJ(1).T = options.tempCushion;
        Wref_INJ(1).sign = 1;
        rng(0)

        dt       = scheduleRef.step.val;
        nstep    = numel(scheduleRef.step.val);
        
        scheduleRef.control = scheduleRef.control(1:scheduleRef.step.control(ninj));
        
        % For the training run, we create a schedule with varying controls
        perturbStep = [ones(1,4), round(.5+(2:(nstep-3))/2)]';
        rng(0)
        scheduleRef = perturbedSimpleSchedule(dt, 'W', Wref, ...           %endret her fra Wref_INJ !!!!
            'pressureFac', pressureFac, 'rateFac', rateFac, 'perturbStep', perturbStep);
        % Pack the simulation problem with the initial state, model, and schedule
        for i = 1:length(scheduleRef.control) 
            scheduleRef.control(i).bc = ref.control(1).bc;
        end
        schedule_inj = scheduleRef;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% production rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nprod > 0
        ind = 1:nprod;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref     = scheduleRef.control.W;

        Wref_PROD = scheduleRef.control.W;
        Wref_PROD(1).type = 'rate';
        Wref_PROD(1).name = 'charge';
        Wref_PROD(1).val = options.rateCharge;
        Wref_PROD(1).T = options.tempCharge;
        Wref_PROD(1).sign = 1;
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

        schedule_bhp = scheduleRef;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% shutting rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nshutt > 0
        ind = 1:nshutt;
        scheduleRef = ref;
        scheduleRef.step.val     = ref.step.val(ind);
        scheduleRef.step.control = ref.step.control(ind);
        Wref     = scheduleRef.control.W;

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


