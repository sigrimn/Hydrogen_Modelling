function problem = createOutputDirectory(name, problem)

% Check and create the directory if it doesn't exist
if ~exist(name, 'dir')
    mkdir(name);
end

% Create subdirectories for wellSols, states, and reports
mkdir(fullfile(name, 'GenericBlackOilModel'));
mkdir(fullfile(name, 'wellSols'));
mkdir(fullfile(name, 'states'));
mkdir(fullfile(name, 'reports'));
% here: simulatePackedProblem with new name and schedule
problem.OutputHandlers.wellSols.dataDirectory = name;
problem.OutputHandlers.states.dataDirectory = name;
problem.OutputHandlers.reports.dataDirectory = name;