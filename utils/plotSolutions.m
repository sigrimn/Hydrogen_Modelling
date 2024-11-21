function varargout = plotSolutions(reference, coarse_solutions, varargin)
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


N = length(reference);

%solRef = struct();

% Assemble solutions
solRef.bhp = zeros(N);
solRef.qOs = zeros(N);
solRef.qGs = zeros(N);

sols = {};

for i = 1:N
    solRef.bhp(i) = reference{i}.bhp;
    solRef.qOs(i) = reference{i}.qOs;
    solRef.qGs(i) = reference{i}.qGs;
end

if length(coarse_solutions) == N
    sols.bhp = zeros(N);
    sols.qOs = zeros(N);
    sols.qGs = zeros(N);
    for j = 1:N
        sols.bhp(j) = coarse_solutions{j}.bhp;
        sols.qOs(j) = coarse_solutions{j}.qOs;
        sols.qGs(j) = coarse_solutions{j}.qGs;
    end
    
else
    for i = 1:length(coarse_solutions)
        sol.bhp = zeros(N);
        sol.qOs = zeros(N);
        sol.qGs = zeros(N);
        for j = 1:N
            sol.bhp(j) = coarse_solutions{i}{j}.bhp;
            sol.qOs(j) = coarse_solutions{i}{j}.qOs;
            sol.qGs(j) = coarse_solutions{i}{j}.qGs;
        end
    
        sols{end+1} = sol;
    end
end
% Generer et sett med farger fra en fargekart
numSolutions = length(coarse_solutions) + 1;
cmap = lines(numSolutions); % Bruker 'lines'-fargekart for ulike farger


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Første subplot: bhp
subplot(3, 1, 1); % 1 rad, 3 kolonner, dette er det første plottet
plot(solRef.bhp, 'Color', cmap(3, :)); % Tildel en unik farge fra cmap 
hold on;
% Plott alle nye løsninger
if length(coarse_solutions) == N
    plot(sols.bhp, 'Color', cmap(1, :));
else
    for i = 1:length(coarse_solutions)
        plot(sols{i}.bhp, 'Color', cmap(i, :)); % Tildel en unik farge fra cmap
    end
end
colormap(cmap);                   % Bruk tilpasset fargekart
colorbar;                         % Fargekart for å vise bergartstyper

% Tilpass tittel og akser for første plott
title('BHP');
xlabel('time t');
ylabel('pressure');
view(2); % 2D-visning
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andre subplott H2 rate
subplot(3, 1, 2); % 1 rad, 3 kolonner, dette er det første plottet
plot(solRef.qOs, 'Color', cmap(3, :)); % Tildel en unik farge fra cmap 
hold on;
% Plott alle nye løsninger
% Plott alle nye løsninger
if length(coarse_solutions) == N
    plot(sols.qOs, 'Color', cmap(1, :));
else
    for i = 1:length(coarse_solutions)
        plot(sols{i}.qOs, 'Color', cmap(i, :)); % Tildel en unik farge fra cmap
    end
end
colorbar;                         % Fargekart for å vise bergartstyper

% Tilpass tittel og akser for andre plott
title('H2 rate');
xlabel('time t');
ylabel('rate');
view(2); % 2D-visning
axis tight;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tredje subplott water rate
subplot(3, 1, 3); % 1 rad, 3 kolonner, dette er det første plottet
plot(solRef.qGs, 'Color', cmap(3, :)); % Tildel en unik farge fra cmap 
hold on;
% Plott alle nye løsninger
if length(coarse_solutions) == N
    plot(sols.qGs, 'Color', cmap(1, :));
else
    for i = 1:length(coarse_solutions)
        plot(sols{i}.qGs, 'Color', cmap(i, :)); % Tildel en unik farge fra cmap
    end
end
    
colormap(cmap);                   % Bruk tilpasset fargekart
colorbar;                         % Fargekart for å vise bergartstyper

% Tilpass tittel og akser for tredje plott
title('H2O rate');
xlabel('time t');
ylabel('pressure');
view(2); % 2D-visning
axis tight;