function varargout = plotWell2d(reference, coarse, reference_schedule, coarse_schedule, varargin)
% Run a schedule for a non-linear physical model using an automatic differention
%
% SYNOPSIS:
%   wellSols = simulateScheduleAD(initState, model, schedule)
%
%   [wellSols, state, report]  = simulateScheduleAD(initState, model, schedule)
%
% DESCRIPTION:
%   This function takes in a valid schedule file (see required parameters)
%   and runs a simulation through all timesteps for a given model and
%   initial state.
%
%   simulateScheduleAD is the outer bookkeeping routine used for running
%   simulations with non-trivial control changes. It relies on the model
%   and (non)linear solver classes to do the heavy lifting.
%
% REQUIRED PARAMETERS:
%       
%
% OPTIONAL PARAMETERS:
%
%   
%                       
% RETURNS:
% plot of reference grid and coarse grid with wells. 
% NOTE:


% Tilpasset fargekart for lysere farger
cmap = [0.6 0.8 1.0;  % Lys blå for type 1
        0.7 1.0 0.7;  % Lys grønn for type 2
        1.0 0.8 0.6]; % Lys oransje for type 3

% Figur for side-ved-side plott
figure;

% Første subplot: modelRef
subplot(1, 2, 1); % 1 rad, 2 kolonner, dette er det første plottet
plotGrid(reference.G, 'FaceColor', 'none'); % Plott rammen av grid først
hold on;
plotCellData(reference.G, reference.rock.regions.saturation); % Fargelegg cellene
colormap(cmap);                   % Bruk tilpasset fargekart
colorbar;                         % Fargekart for å vise bergartstyper

% Legg til brønnen i modelRef
if ~isempty(reference_schedule.control(1).W)
    plotGrid(reference.G, reference_schedule.control(1).W.cells, 'EdgeAlpha', 0.8);
end

% Tilpass tittel og akser for første plott
title('Rock Types and Well (modelRef)');
xlabel('x [m]');
ylabel('y [m]');
view(2); % 2D-visning
axis tight;

% Andre subplot: modelCoarse
subplot(1, 2, 2); % 1 rad, 2 kolonner, dette er det andre plottet
plotGrid(coarse.G, 'FaceColor', 'none'); % Plott rammen av grid først
hold on;
plotCellData(coarse.G, coarse.rock.regions.saturation); % Fargelegg cellene
colormap(cmap);                   % Bruk tilpasset fargekart
colorbar;                         % Fargekart for å vise bergartstyper

% Legg til brønnen i modelCoarse
if ~isempty(coarse_schedule.control(1).W)
    plotGrid(coarse.G, coarse_schedule.control(1).W.cells, 'EdgeAlpha', 0.8);
end

% Tilpass tittel og akser for andre plott
title('Rock Types and Well (modelCoarse)');
xlabel('x [m]');
ylabel('y [m]');
view(2); % 2D-visning
axis tight;
