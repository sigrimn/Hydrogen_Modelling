function varargout = fullSimulationComparison(wsRef, wsCoarse1, wsCoarse2, varargin)

%creating vectors to contain the data
N = length(wsRef);
reference_bhp = zeros(N);
reference_qOs = zeros(N);
reference_qGs = zeros(N);

new_bhp1 = zeros(N);
new_qOs1 = zeros(N);
new_qGs1 = zeros(N);

new_bhp2 = zeros(N);
new_qOs2 = zeros(N);
new_qGs2 = zeros(N);


for i = 1:2518
    reference_bhp(i) = wsRef{i}.bhp;
    reference_qOs(i) = wsRef{i}.qOs;
    reference_qGs(i) = wsRef{i}.qGs;

    new_bhp1(i) = wsCoarse1{i}.bhp;
    new_qOs1(i) = wsCoarse1{i}.qOs;
    new_qGs1(i) = wsCoarse1{i}.qGs;

    new_bhp2(i) = wsCoarse2{i}.bhp;
    new_qOs2(i) = wsCoarse2{i}.qOs;
    new_qGs2(i) = wsCoarse2{i}.qGs;
end


% Plotting
%time = 1:N;  % Create a time vector corresponding to the length of A and B

% Plotting reference and fitted full models
figure; % Create a new figure
plot(reference_bhp, '-b');  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_bhp1, '-r');  % Plot B with a red line
hold on;
plot(new_bhp2, '-g');  % Plot B with a red line
% Adding legends
legend('Reference model', 'Fitted model BFGS','Fitted model, L-M');

% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('barsa');  % Label for the y-axis
title('Comparison of reference model and fitted model bottom hole pressure');  % Title for the plot



% Optional: Grid for better visualization
grid on;
hold off;
%%
% Plotting reference and fitted full models
figure; % Create a new figure
plot(reference_qOs, '-b', 'LineWidth', 1.5);  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_qOs1, '-r', 'LineWidth', 1.5);  % Plot B with a red line
hold on;
plot(new_qOs2, '-g');  % Plot B with a red line

% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('m^3/day');  % Label for the y-axis
title('Comparison of reference model and fitteed model water flow');  % Title for the plot

% Adding legends
legend('Reference model', 'Fitted model BFGS', 'Fitted model, L-M');

% Optional: Grid for better visualization
grid on;
hold off;

%%
% Plotting reference and fitted full models
figure; % Create a new figure

plot(reference_qGs, '-b', 'LineWidth', 1.5);  % Plot A with a blue line
hold on;  % Hold the current plot to overlay the next plot
plot(new_qGs1, '-r', 'LineWidth', 1.5);  % Plot B with a red line
hold on;  % Hold the current plot to overlay the next plot
plot(new_qGs2, '-r', 'LineWidth', 1.5);  % Plot B with a red line

% Adding labels and title
xlabel('time');  % Label for the x-axis
ylabel('m^3/day');  % Label for the y-axis
title('Comparison of reference model and fitteed model hydrogen flow');  % Title for the plot

legend('Reference model', 'Fitted model BFGS','Fitted model, L-M');


%legend('Reference model', 'Fitted model, L-M', 'b', 'r');
%ah1 = axes('position',get(gca,'position'),'visible','off');
%legend(ah1, [h1(2) h2(2)], {'Test1','Test2'}, 'Reference model', 'Fitted model, L-M');

% Optional: Grid for better visualization
grid on;
hold off;