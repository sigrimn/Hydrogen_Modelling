clear all;
close all;

%Open variables from 'data.mat'
load('data_121124.mat');
% Now, a and b are available in the workspace


% Antall parametere
numParamsFluid = 2000; % Antall "FluidVolume" parametere
numParamsTrans = 7000; % Antall "Transmissibilities" parametere
numParams = numParamsFluid + numParamsTrans;

% Skalerte verdier for initial og final X
initialX = rand(1, numParams) * 0.001; % Tilfeldige initialverdier
finalX = initialX + (rand(1, numParams) - 0.5) * 0.0002; % Endelige verdier etter optimalisering

% Lower bound (antatt konstant verdi for eksempel)
lowerBound = zeros(1, numParams); % Kan være justert etter behov

% X-aksen for hver parametergruppe
xAxis = 1:numParams;
fluidEnd = numParamsFluid; % Siste indeks for "FluidVolume"
transEnd = numParams; % Siste indeks for "Transmissibilities"

% Opprett plottet
figure;
hold on;
scatter(xAxis, initialX, 10, 'b', 'filled'); % Initial X verdier (blå)
scatter(xAxis, finalX, 10, 'orange', 'filled'); % Final X verdier (oransje)
plot(xAxis, lowerBound, 'cyan', 'LineWidth', 1.5); % Lower bound (cyan linje)

% Legg til tekstetiketter for parametergruppene
text(fluidEnd/2, -0.00005, 'FluidVolume', 'HorizontalAlignment', 'center', 'FontSize', 12);
text((fluidEnd + transEnd)/2, -0.00005, 'Transmissibilities', 'HorizontalAlignment', 'center', 'FontSize', 12);

% Formatering av plottet
title('Scaled parameters');
ylabel('Scaled value');
xlabel('');
legend({'Initial X', 'Final X', 'Lower bound'}, 'Location', 'northeast');
grid on;
hold off;
